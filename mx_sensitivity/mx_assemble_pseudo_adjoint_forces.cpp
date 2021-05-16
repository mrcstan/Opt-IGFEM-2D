/*
Created by Marcus Tan on 12/23/2014
Modified on 2/21/2015
Copyright 2014 University of Illinois 
Purpose: this function calculates the pseudo and adjoint forces
         and the remainder term in the gradient of the objective
         function wrt design parameters. 
         see assemble_pseudo_adjoint_forces for formulation
NOTE:
    i): For nonlinear radiation analysis, temperature unit must be Kelvin
Input:
    (1): objOpt: structural variable containing information for type of sensitivity analysis
                 see sensitivy.h for members of this variable
    (2): nodeCoords:
    (3): elemNodes:
    (4): elemHeatSource:
    (5): gauss1:
    (6): parent:
    (7): channels: channels.D2mcf: nDesignParams x nChannels matrix,
                   each entry corresponding to the derivative of the
                   heat capacity x mass flow rate vs the corresponding
                   design parameter
	(8): Neumann:
    (9): UUR:
    (10): supg: true or false
Output:
    (1): Fpseudo: a nDesignParams x nNodes pseudo force matrix 
                 (also called trans(Pps) in notes)
                 for column-wise write access, I defined it that way 
                 instead of the transpose
    (2): Fadj: an adjoint force vector of size nNodes 
               (also called trans(D1_PI) in notes)
    (3): gradObjRem: a vector of size nDesignParams of the remainder 
                     term of the gradient of the objective function 
                     wrt design parameters (also called D3_PI in notes)
    (4): objVal: objective function values
    (5): Fadj1norm: the adjoint force vector for 1-norm temperature
    (6): grad1normRem: remainder term of the gradient of the 1-norm temperature
    (7): oneNormVal: 1-norm temperature value
*/
#include <iostream>
#include <cstddef> // for size_t
#include <math.h>
#include <algorithm>
#include "sensitivity.h"
#include "armaMex.hpp"

using namespace igfem;

enum Inputs {OBJ_OPTIONS,NODE_COORDS,ELEM_NODES,ELEM_HEAT,CONVECTION,
             GAUSS,PARENT,CHANNEL_NETWORK,NEUMANN,UNREDUCED_SOLN,
             SUPG,CALC_GRAD,IN_LAST};
const int nInputs = IN_LAST - OBJ_OPTIONS;


enum ChFields {ISTRIANGLE,LOC_PA_NODES,LOC_PA_EN_NODES,LOC_EN_NODES,
               CHANNEL_NUM,CHANNEL_NODES,CHANNEL_LOC_NODES,CH_LAST};
const int nChFields = CH_LAST - ISTRIANGLE;

enum PaFields {PA_TYPE,PA_NODES,PA_LOC_CHAN_NODES,PA_CHAN_NUM,PA_CONDUCTIVITY,
               PA_CHILD,PA_VEL,PA_DESIGN_PARAM_NUM,PA_LAST};
const int nPaFields = PA_LAST - PA_TYPE;

enum Outputs {PSEUDO_FORCE,ADJOINT_FORCE,GRAD_OBJ_REMAINDER,OBJ_VAL,
              ADJOINT_FORCE_1NORM, GRAD_1NORM_REMAINDER,ONENORM,OUT_LAST};
const int nOutputs = OUT_LAST - PSEUDO_FORCE;


void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
    /* check for proper number of arguments */
    if (nrhs < nInputs) 
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:nrhs",
                          "%i inputs required.",nInputs);
    
    if (nlhs > nOutputs) 
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:nlhs",
                          "At most %i outputs can be requested.",nOutputs);
     
    if (!mxIsStruct(prhs[OBJ_OPTIONS]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:OBJ_OPTIONS",
                          "Input %i, objective options must be structural",OBJ_OPTIONS+1);
    
    if (!mxIsDouble(prhs[NODE_COORDS]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:NODE_COORDS",
                          "Input %i, nodeCoords must be double",NODE_COORDS+1);
    
    if (!mxIsDouble(prhs[ELEM_NODES]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:ELEM_NODES",
                          "Input %i, elemNodes must be double",ELEM_NODES+1);

    if (!mxIsDouble(prhs[ELEM_HEAT]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:ELEM_HEAT",
                          "Input %i, elemHeatSource must be double",ELEM_HEAT+1);
    
    if (!mxIsStruct(prhs[CONVECTION]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:CONVECTION",
                          "Input %i, convection must be structural", CONVECTION);
    
    if (!mxIsStruct(prhs[GAUSS]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:GAUSS",
                          "Input %i, gauss must be structural",GAUSS+1);
    
    if (!mxIsStruct(prhs[PARENT]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:PARENT",
                          "Input %i, parent must be structural",PARENT+1);
    
    if (!mxIsStruct(prhs[CHANNEL_NETWORK]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:CHANNEL_NETWORK",
                          "Input %i, channel network must be structural",CHANNEL_NETWORK+1);
   
   if (!mxIsStruct(prhs[NEUMANN]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:NEUMANN",
                          "Input %i, Neumann must be structural",NEUMANN+1);
    
 
    if (!mxIsDouble(prhs[UNREDUCED_SOLN]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:UNREDUCED_SOLN",
                          "Input %i, UUR must be double",UNREDUCED_SOLN+1);

    if (!mxIsLogicalScalar(prhs[SUPG]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:SUPG",
                          "Input %i, supg must be scalar boolean",SUPG+1);
    
    if (!mxIsLogicalScalar(prhs[CALC_GRAD]))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:CALC_GRAD",
                          "Input %i, calcGrad must be scalar boolean",CALC_GRAD+1);
                       
    bool supg = mxIsLogicalScalarTrue(prhs[SUPG]); 

    bool calcGrad = mxIsLogicalScalarTrue(prhs[CALC_GRAD]);
	
    mxArray* fieldPtr;
    int fieldNum;

    ObjOptions objOpt; 
    fieldNum = mxGetFieldNumber(prhs[OBJ_OPTIONS],"normp");
    if (fieldNum < 0)
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_obj_normp",
                          "objective function pnorm value not found");
    fieldPtr = mxGetFieldByNumber(prhs[OBJ_OPTIONS],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:obj_normp_empty",
                         "objective function pnorm value not provided");
    if(mxIsDouble(fieldPtr))
    {
        objOpt.normp = mxGetScalar(fieldPtr);
        if (objOpt.normp <= 0)
            mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:obj_normp_nonpositive",
                              "objective function pnorm value must be positive");
    }
    else
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:obj_normp_not_double",
                          "objective function pnorm value must be of type double");
    
    fieldNum = mxGetFieldNumber(prhs[OBJ_OPTIONS],"offset");
    if (fieldNum < 0)
        objOpt.offset = 0.0;
    else
    {
        fieldPtr = mxGetFieldByNumber(prhs[OBJ_OPTIONS],0,fieldNum);
        if(!fieldPtr || mxIsEmpty(fieldPtr))
            objOpt.offset = 0.0;
        else if (mxIsDouble(fieldPtr))
            objOpt.offset = mxGetScalar(fieldPtr);
        else
            objOpt.offset = 0.0;
    }
    //mexPrintf("offset = %g\n",objOpt.offset);

    fieldNum = mxGetFieldNumber(prhs[OBJ_OPTIONS],"nDesignParams");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_nDesignParams",
                           "objective function number of design parameters not provided");
    fieldPtr = mxGetFieldByNumber(prhs[OBJ_OPTIONS],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:nDesignParams_empty",
                            "objective function number of design parameters is empty");
    objOpt.nDesignParams = mxGetScalar(fieldPtr);
    
    fieldNum = mxGetFieldNumber(prhs[OBJ_OPTIONS],"calcOneNorm");
    if (fieldNum < 0)
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_obj_calcOneNorm",
                          "objective function calcOneNorm boolean not found");
    fieldPtr = mxGetFieldByNumber(prhs[OBJ_OPTIONS],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:obj_calcOneNorm_empty",
                         "objective function calcOneNorm boolean not provided");
    if(mxIsLogical(fieldPtr))
    {
        objOpt.calcOneNorm = mxIsLogicalScalarTrue(fieldPtr);
    }
    
    int intDomain;
    fieldNum = mxGetFieldNumber(prhs[OBJ_OPTIONS],"intDomainType");
    if (fieldNum < 0)
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_int_domain_type",
                          "objective function intDomainType not found");
    fieldPtr = mxGetFieldByNumber(prhs[OBJ_OPTIONS],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:obj_int_domain_empty",
                         "objective function intDomainType is empty");
    if(mxIsDouble(fieldPtr))
    {
        intDomain = mxGetScalar(fieldPtr);
        if (intDomain == 0)
            objOpt.intDomain = WHOLE;
        else if (intDomain == 1)
            objOpt.intDomain = CHANNEL;
        else
            mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:obj_int_domain_type_unknown",
                              "unknown objective function intDomainType");
    }
    else
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:obj_int_domain_type_not_double",
                          "objective function intDomainType must be of type double");

    //mexPrintf("get nodal coordinates \n");
    if (mxGetM(prhs[NODE_COORDS]) != 2)
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:nodeCoordsRow",
                         "number of rows of nodeCoords must be 2");
    arma::mat nodeCoords = armaGetPr(prhs[NODE_COORDS]);

    //mexPrintf("get elem nodes \n"):
    if (mxGetM(prhs[ELEM_NODES]) != N_ORIGINAL_NODES)
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:elemNodesRow",
                         "number of rows of elemNodes must be %i",N_ORIGINAL_NODES);
    arma::umat elemNodes = arma::conv_to<arma::umat>::from(armaGetPr(prhs[ELEM_NODES])-1);

    //mexPrintf("get element heat source \n");
    arma::vec elemHeatSource;
    if (mxGetM(prhs[ELEM_HEAT]) == 1)
        elemHeatSource = arma::trans(armaGetPr(prhs[ELEM_HEAT]));  
    else if (mxGetN(prhs[ELEM_HEAT]) == 1)
        elemHeatSource = armaGetPr(prhs[ELEM_HEAT]);
    else 
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:elemHeatSource_vector",
                          "elemHeatSource must be row or column vector");        
	
    // validating and getting gauss integration points
    gauss gauss1;
    //mexPrintf("validating and getting gauss integration points \n");
    // element integration
    fieldNum = mxGetFieldNumber(prhs[GAUSS],"elem");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_gauss_elem",
                           "gauss points for element integration not found");
    fieldPtr = mxGetFieldByNumber(prhs[GAUSS],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:gauss_elem_empty",
                            "gauss points for element integration must be > 0");
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:gauss_elem_not_double",
                            "gauss points for element integration must be of type double");
    gauss1.elem = armaGetPr(fieldPtr);  
    if(gauss1.elem.n_rows != nodeCoords.n_rows+1)
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:gauss_elem_wrong_dim",
                            "wrong dimension for element integration");
    
    // quadrilateral element integration
    fieldNum = mxGetFieldNumber(prhs[GAUSS],"quadElem");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_gauss_quad_elem",
                           "gauss points for quad element integration not found");
    fieldPtr = mxGetFieldByNumber(prhs[GAUSS],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:gauss_quad_elem_empty",
                            "gauss points for quad element integration must be > 0");
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:gauss_quad_elem_not_double",
                            "gauss points for quad element integration must be of type double");
    gauss1.quadElem = armaGetPr(fieldPtr);  
    if(gauss1.quadElem.n_rows != nodeCoords.n_rows+1)
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:gauss_quad_elem_wrong_dim",
                            "wrong dimension for quad element integration");

     // line integration
    fieldNum = mxGetFieldNumber(prhs[GAUSS],"line");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_gauss_line",
                           "gauss points for line integration not found");
    fieldPtr = mxGetFieldByNumber(prhs[GAUSS],0,fieldNum);
    if (!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:gauss_line_empty",
                            "gauss points for line integration must be > 0");
    
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:gauss_line_not_double",
                            "gauss points for line integration must be of type double");
   
    gauss1.line = armaGetPr(fieldPtr);
    if (gauss1.line.n_rows != 2)
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:gauss_line_wrong_dim",
                            "wrong dimension for line integration");
							
    // validating and getting parent data
    //mexPrintf("validating and getting parent data \n");   
    std::size_t nElems = mxGetNumberOfElements(prhs[PARENT]);
    arma::field<parent> parents(nElems);
    
    int paFieldNum[nPaFields];
    for (int i = 0; i < nPaFields; i++)
        paFieldNum[i] = -1;
    
    mxArray* chFieldPtr;
    int chFieldNum[nChFields]; 
    for (int i = 0; i < nChFields; i++)
        chFieldNum[i] = -1;

    // parent type
    paFieldNum[PA_TYPE] = mxGetFieldNumber(prhs[PARENT],"type");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_type",
                           "parent type field not found");
    // parent local channel nodes
    paFieldNum[PA_LOC_CHAN_NODES] = mxGetFieldNumber(prhs[PARENT],"channelLocNodes");
    if (fieldNum < 0)
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_channelLocNodes",
                          "parent(i).channelLocNodes field not found");
    // parent channel numbers
    paFieldNum[PA_CHAN_NUM] = mxGetFieldNumber(prhs[PARENT],"channelNum");
    if (fieldNum < 0)
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_channelNum",
                          "parent(i).channelNum field not found");
    
    bool getPaChFieldNum = true;
    int parentType = 0;
    for (std::size_t j = 0; j < parents.n_elem; j++)
    {
        fieldPtr = mxGetFieldByNumber(prhs[PARENT],j,paFieldNum[PA_TYPE]);        
        if (!fieldPtr || mxIsEmpty(fieldPtr)) 
        {
            mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:parent_type_not_found",
                              "parent_type_not_found for element %i\n",j);
        }
        parentType = mxGetScalar(fieldPtr);
        if (parentType == 1) // regular element with channels    
        {
            fieldPtr = mxGetFieldByNumber(prhs[PARENT],j,paFieldNum[PA_LOC_CHAN_NODES]);
            if (fieldPtr && !mxIsEmpty(fieldPtr))
            {
                parents(j).channelLocNodes = arma::conv_to<arma::umat>::from(armaGetPr(fieldPtr)-1);
                if(parents(j).channelLocNodes.n_rows != 2)
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:parent_channelLocNodes_no_2_rows",
                                      "parent(i).channelLocNodes must be a 2 x number of channelNodes matrix");
            }
            fieldPtr = mxGetFieldByNumber(prhs[PARENT],j,paFieldNum[PA_CHAN_NUM]);
            if (fieldPtr && !mxIsEmpty(fieldPtr))
            {
                parents(j).channelNum = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
            }

        }
        else if (parentType == 2) // IGFEM element
        {
            //mexPrintf("parent %i \n",j);
            parents(j).type = IGFEM;
            if (getPaChFieldNum)
            {
                paFieldNum[PA_NODES] = mxGetFieldNumber(prhs[PARENT],"nodes");
                if (paFieldNum[PA_NODES] < 0)
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_node",
                                      "parent node field not found");
                
                paFieldNum[PA_CONDUCTIVITY] = mxGetFieldNumber(prhs[PARENT],"conductivity");
                if (paFieldNum[PA_CONDUCTIVITY] < 0)
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_conductivity",
                                      "parent(i).conductivity field not found");
   
                paFieldNum[PA_CHILD] = mxGetFieldNumber(prhs[PARENT],"child");
                if (paFieldNum[PA_CHILD] < 0)
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_child",
                                      "parent(i).child field not found");
             
                if (calcGrad)
                {
                    paFieldNum[PA_VEL] = mxGetFieldNumber(prhs[PARENT],"vel");
                    if (paFieldNum[PA_VEL] < 0)
                        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_vel",
                                          "parent(i).vel field not found");
                    
                    paFieldNum[PA_DESIGN_PARAM_NUM] = mxGetFieldNumber(prhs[PARENT],"designParamNum");
                    if (paFieldNum[PA_DESIGN_PARAM_NUM] < 0)
                        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_designParamNum",
                                          "parent(i).designParamNum field not found");
                }
                // get child field numbers
                fieldPtr = mxGetFieldByNumber(prhs[PARENT],j,paFieldNum[PA_CHILD]);  
                if (!fieldPtr || !mxIsStruct(fieldPtr))
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:parent_child_not_struct",
                                    "parent.child is not initialized or not of type struct");

                // isTriangle
                chFieldNum[ISTRIANGLE] = mxGetFieldNumber(fieldPtr,"isTriangle");
                if (chFieldNum[ISTRIANGLE] < 0)
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_child_isTriangle",
                                      "parent.child.isTriangle field not found");

                // locPaNodes
                chFieldNum[LOC_PA_NODES] = mxGetFieldNumber(fieldPtr,"locPaNodes");
                if (chFieldNum[LOC_PA_NODES] < 0)
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_child_locPaNodes",
                                      "parent.child.locPaNodes field not found");
       
                // locPaEnNodes
                chFieldNum[LOC_PA_EN_NODES] = mxGetFieldNumber(fieldPtr,"locPaEnNodes");
                if (chFieldNum[LOC_PA_EN_NODES] < 0)
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_child_locPaEnNodes",
                                       "parent.child.locPaEnNodes field not found"); 

                // locEnNodes
                chFieldNum[LOC_EN_NODES] = mxGetFieldNumber(fieldPtr,"locEnNodes");
                if (chFieldNum[LOC_EN_NODES] < 0)
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_child_locEnNodes",
                                       "parent.child.locEnNodes field not found");

				// channelNum
				chFieldNum[CHANNEL_NUM] = mxGetFieldNumber(fieldPtr,"channelNum");
				if (chFieldNum[CHANNEL_NUM] < 0)
					mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_child_channelNum",
									   "parent.child.channelNum field not found");


				// channelNodes
				chFieldNum[CHANNEL_NODES] = mxGetFieldNumber(fieldPtr,"channelNodes");
				if (chFieldNum[CHANNEL_NODES] < 0)
					mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_child_channelNodes",
									   "parent.child.channelNodes field not found");
				
				chFieldNum[CHANNEL_LOC_NODES] = mxGetFieldNumber(fieldPtr,"channelLocNodes");
				if (chFieldNum[CHANNEL_LOC_NODES] < 0)
					mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_parent_child_channelLocNodes",
									   "parent.child.channelLocNodes field not found");	
                getPaChFieldNum = false;
            }
            
            // parent nodes
            fieldPtr = mxGetFieldByNumber(prhs[PARENT],j,paFieldNum[PA_NODES]);
            if (fieldPtr && !mxIsEmpty(fieldPtr))
                parents(j).nodes = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
            

            // conductivity
            fieldPtr = mxGetFieldByNumber(prhs[PARENT],j,paFieldNum[PA_CONDUCTIVITY]);
            if (fieldPtr && !mxIsEmpty(fieldPtr))
                parents(j).conductivity = mxGetScalar(fieldPtr);
            
            if (calcGrad)
            {
                // velocity
                fieldPtr = mxGetFieldByNumber(prhs[PARENT],j,paFieldNum[PA_VEL]);
                if (fieldPtr && !mxIsEmpty(fieldPtr))
                {
                    if (mxGetNumberOfDimensions(fieldPtr) == 3)
                    {
                        parents(j).vel = armaGetCubePr(fieldPtr);
                    }
                    else if (mxGetNumberOfDimensions(fieldPtr) == 2)
                    {
                        parents(j).vel.set_size(mxGetM(fieldPtr),mxGetN(fieldPtr),1);
                        parents(j).vel.slice(0) = armaGetPr(fieldPtr);
                    }
                    else
                        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:parent_vel_n_dims",
                                           "parent.vel must have either 2 or 3 dimensions");
                       
                }
                // design parameter numbers
                fieldPtr = mxGetFieldByNumber(prhs[PARENT],j,paFieldNum[PA_DESIGN_PARAM_NUM]);
                if (fieldPtr && !mxIsEmpty(fieldPtr))
                    parents(j).designParamNum = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
            }

            // children
            std::size_t nChilds = 0;
            fieldPtr = mxGetFieldByNumber(prhs[PARENT],j,paFieldNum[PA_CHILD]);  
            if (fieldPtr && mxIsStruct(fieldPtr) && !mxIsEmpty(fieldPtr))
            {
                nChilds = mxGetNumberOfElements(fieldPtr);
                parents(j).children.set_size(nChilds);
            }
            for (std::size_t k = 0; k < nChilds; k++)
            {
                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[ISTRIANGLE]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                {
                    if (mxIsLogicalScalarTrue(chFieldPtr))
                        parents(j).children(k).shape = TRIANGLE;
                    else
                        parents(j).children(k).shape = QUADRILATERAL;
                }
                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[LOC_PA_NODES]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                {
                    parents(j).children(k).locPaNodes = arma::conv_to<arma::uvec>::from(armaGetPr(chFieldPtr)-1);
                }
                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[LOC_PA_EN_NODES]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                    parents(j).children(k).locPaEnNodes = arma::conv_to<arma::uvec>::from(armaGetPr(chFieldPtr)-1);

                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[LOC_EN_NODES]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                    parents(j).children(k).locEnNodes = arma::conv_to<arma::uvec>::from(armaGetPr(chFieldPtr)-1);   

			   	 chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[CHANNEL_NUM]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                    parents(j).children(k).channelNum = arma::conv_to<arma::uvec>::from(armaGetPr(chFieldPtr)-1);                

                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[CHANNEL_NODES]);
                             
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                {
                    std::size_t nChannels = mxGetNumberOfElements(chFieldPtr);
                    mxArray *cellElemPtr;
                    parents(j).children(k).channelNodes.set_size(2,nChannels);
                    arma::uvec channelNodes;
                    for (std::size_t i = 0; i < nChannels; i++)
                    {
                        cellElemPtr = mxGetCell(chFieldPtr,i);
                        if (cellElemPtr && !mxIsEmpty(cellElemPtr))
                        {
                            channelNodes = arma::conv_to<arma::uvec>::from(armaGetPr(cellElemPtr)-1);
                            parents(j).children(k).channelNodes(0,i) = channelNodes(0);
                            parents(j).children(k).channelNodes(1,i) = channelNodes(channelNodes.n_elem-1);
                        }
                    }
                }

                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[CHANNEL_LOC_NODES]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                {
                    parents(j).children(k).channelLocNodes 
                        = arma::conv_to<arma::umat>::from(armaGetPr(chFieldPtr)-1);
                }			

            }


        } // parent.type == 2
    }
    
	// validating channel data
    //mexPrintf("validating and getting channel data \n");
    fieldNum = mxGetFieldNumber(prhs[CHANNEL_NETWORK],"mcf");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_channel_mcf",
                           "channel mcf field not found");
    fieldPtr = mxGetFieldByNumber(prhs[CHANNEL_NETWORK],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:channel_mcf_empty",
                             "channel mcf field empty");   
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:channel_mcf_not_double",
                            "channel must be of type double"); 
							
	chanNetwork channels;
    channels.mcf = armaGetPr(fieldPtr); 
    
    fieldNum = mxGetFieldNumber(prhs[CHANNEL_NETWORK],"model");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_channel_model",
                           "channel model field not found");
    fieldPtr = mxGetFieldByNumber(prhs[CHANNEL_NETWORK],0,fieldNum);  
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:channel_model_empty",
                             "channel model field empty");   
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:channel_model_not_double",
                            "channel model must be of type double"); 

    
    int modelNum = mxGetScalar(fieldPtr); 
    if (modelNum  == 1)
        channels.model = MEAN_TEMP;
    else if (modelNum == 2)
    {
         mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:const_heat_not_yet",
                            "constant heat flux model not yet implemented"); 
    }
    else
		 mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:unknown_mode",
                            "unknown channel model"); 
    
    fieldNum = mxGetFieldNumber(prhs[CHANNEL_NETWORK],"D2mcf");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:no_channel_D2mcf",
                           "channel D2mcf field not found");
    fieldPtr = mxGetFieldByNumber(prhs[CHANNEL_NETWORK],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:channel_D2mcf_empty",
                             "channel D2mcf field empty");   
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:channel_mcf_not_double",
                            "channel must be of type double"); 
    //channels.D2mcf = armaGetSparseMatrix(fieldPtr,true); // sort locations
    channels.D2mcf = arma::trans(armaGetPr(fieldPtr)); 
    //channels.D2mcf.print("D2mcf = ");
	
	mexPrintf("validating and getting Neumann info \n");  
    arma::field<Neumann> elemNeumann(elemNodes.n_cols);

    arma::uvec neumannElems,neumannSurfs;
    arma::vec neumannVals;

    fieldNum = mxGetFieldNumber(prhs[NEUMANN],"heatFlux_elem");
    if (fieldNum > 0)
    {
        // heatFlux_elem
        fieldPtr = mxGetFieldByNumber(prhs[NEUMANN],0,fieldNum);        
        if (fieldPtr && !mxIsEmpty(fieldPtr))
        {
            neumannElems = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
            
            
            
            // heatFlux_surface
            fieldNum = mxGetFieldNumber(prhs[NEUMANN],"heatFlux_surface");
            if (fieldNum < 0)
                mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:heatFlux_surface",
                                   "neumann heatFlux_surface field not found");
            else
            {
                fieldPtr = mxGetFieldByNumber(prhs[NEUMANN],0,fieldNum);
                if (fieldPtr && !mxIsEmpty(fieldPtr))
                    neumannSurfs = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
                else
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:heatFlux_surface",
                                      "neumann heatFlux_surface field not found");
            }

            // heatFlux_value
            fieldNum = mxGetFieldNumber(prhs[NEUMANN],"heatFlux_value");
            if (fieldNum < 0)
                mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:heatFlux_value",
                                   "neumann heatFlux_value field not found");
            else
            {
                fieldPtr = mxGetFieldByNumber(prhs[NEUMANN],0,fieldNum);
                if (fieldPtr && !mxIsEmpty(fieldPtr))
                    neumannVals = armaGetPr(fieldPtr);
                else
                    mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:heatFlux_value",
                                      "neumann heatFlux_value field not found");
            }
        }
        for (std::size_t i=0; i < neumannElems.n_elem; i++)
        {
            elemNeumann(neumannElems(i)).surf = neumannSurfs(i);
            elemNeumann(neumannElems(i)).val = neumannVals(i);
        }
    }
	
    //mexPrintf("getting unreduced solution\n");
    arma::vec UUR;
    if (mxGetM(prhs[UNREDUCED_SOLN]) == 1)
        UUR = arma::trans(armaGetPr(prhs[UNREDUCED_SOLN]));  
    else if (mxGetN(prhs[UNREDUCED_SOLN]) == 1)
        UUR = armaGetPr(prhs[UNREDUCED_SOLN]);
    else 
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:uur_vector",
                          "UUR must be row or column vector");        
	
	// the number of elements of uur must be the same as the number of nodes
    if (nodeCoords.n_cols > UUR.n_elem)
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:uur_n_elem_no_of_nodes",
                          "UUR.n_elem must be greater than the number of nodes");

    
    
    Convection convect;
    fieldPtr = mxGetField(prhs[CONVECTION],0,"coef");
    if (fieldPtr && !mxIsEmpty(fieldPtr) && mxIsDouble(fieldPtr))
        convect.coef = mxGetScalar(fieldPtr);
    else
    {
        mexWarnMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:convection_coef",
                            "convection coef not found");
        convect.coef = 0.0;
    }
    fieldPtr = mxGetField(prhs[CONVECTION],0,"epsSB");
    if (fieldPtr && !mxIsEmpty(fieldPtr) && mxIsDouble(fieldPtr))
        convect.epsSB = mxGetScalar(fieldPtr);
    else
    {
        mexWarnMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:epsSB",
                           "emissivity*Stefan-Boltzmann not found");
        convect.epsSB = 0.0;
    }
    //mexPrintf("epsSB = %g\n",convect.epsSB);

    fieldPtr = mxGetField(prhs[CONVECTION],0,"Tamb");
    if (fieldPtr && !mxIsEmpty(fieldPtr) && mxIsDouble(fieldPtr))
        convect.Tamb = mxGetScalar(fieldPtr);
    else
    {
        mexWarnMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:convection_Tamb",
                           "convection ambient temperature not found");
        convect.Tamb = 0.0;
    }
    
    fieldPtr = mxGetField(prhs[CONVECTION],0,"linearRad");
    if (fieldPtr && !mxIsEmpty(fieldPtr))
        convect.linearRad = mxIsLogicalScalarTrue(fieldPtr);
    else
    {
        mexWarnMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:convection_linearize_radiation",
                           "linearize radiation flag not found");
        convect.linearRad = true;
    }
    
    fieldPtr = mxGetField(prhs[CONVECTION],0,"TunitNum");
    int TunitNum = mxGetScalar(fieldPtr);
    if (TunitNum  == KELVIN)
    {
        convect.Tunit = KELVIN;
        if (convect.linearRad)
            convect.coef += 4.0*convect.epsSB*pow(convect.Tamb,3.0);
    }
    else if (TunitNum == CELSIUS)
    {
        convect.Tunit = CELSIUS;
        if (convect.linearRad)
            convect.coef += 4.0*convect.epsSB*pow(convect.Tamb+273.15,3.0);
    }
    else
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:convection_TunitNum",
                           "convection temperature unit number not recognized");

    /*
    arma::vec UP;
    if (mxGetM(prhs[PRESCRIBED_TEMP]) == 1)
        UP = arma::trans(armaGetPr(prhs[PRESCRIBED_TEMP]));  
    else if (mxGetN(prhs[PRESCRIBED_TEMP]) == 1)
        UP = armaGetPr(prhs[PRESCRIBED_TEMP]);
    else 
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:UP_vector",
                          "UP must be row or column vector");        
	
	// the number of elements of uur must be the same as the number of nodes
    if (nodeCoords.n_cols > UP.n_elem)
        mexErrMsgIdAndTxt("mx_assemble_pseudo_adjoint_forces:UP_n_elem_no_of_nodes",
                          "UP.n_elem must be greater than the number of nodes");
    */

    arma::mat Fpseudo;
    arma::vec Fadj, gradObjRem;
    double objVal; 
    arma::vec Fadj1norm, grad1normRem;
    double oneNormVal; 
    
    mexPrintf("assembling pseudo adjoint forces\n");
    assemble_pseudo_adjoint_forces(Fpseudo,
                                   Fadj,
                                   gradObjRem,
                                   objVal,
                                   Fadj1norm,
                                   grad1normRem,
                                   oneNormVal,
                                   objOpt,
                                   nodeCoords,
                                   elemNodes,
                                   elemHeatSource,
                                   convect,
                                   gauss1,
                                   parents,
                                   channels,
                                   elemNeumann,
                                   UUR,
                                   supg,
                                   calcGrad);
    
    plhs[PSEUDO_FORCE] = armaCreateMxMatrix(Fpseudo.n_rows,Fpseudo.n_cols);
    plhs[ADJOINT_FORCE] = armaCreateMxMatrix(Fadj.n_elem,1);
    plhs[GRAD_OBJ_REMAINDER] = armaCreateMxMatrix(gradObjRem.n_elem,1);
    armaSetPr(plhs[PSEUDO_FORCE],Fpseudo);
    armaSetPr(plhs[ADJOINT_FORCE],Fadj);
    armaSetPr(plhs[GRAD_OBJ_REMAINDER],gradObjRem);
    plhs[OBJ_VAL] = mxCreateDoubleScalar(objVal);
    
    if (objOpt.calcOneNorm) 
    {
        plhs[ADJOINT_FORCE_1NORM] = armaCreateMxMatrix(Fadj1norm.n_elem,1);
        plhs[GRAD_1NORM_REMAINDER] = armaCreateMxMatrix(grad1normRem.n_elem,1);
        armaSetPr(plhs[ADJOINT_FORCE_1NORM],Fadj1norm);
        armaSetPr(plhs[GRAD_1NORM_REMAINDER],grad1normRem);
        plhs[ONENORM] = mxCreateDoubleScalar(oneNormVal);
    }
    else
    {
        plhs[ADJOINT_FORCE_1NORM] = mxCreateDoubleScalar(0);
        plhs[GRAD_1NORM_REMAINDER] = mxCreateDoubleScalar(0);
        plhs[ONENORM] = mxCreateDoubleScalar(0);
    }
}

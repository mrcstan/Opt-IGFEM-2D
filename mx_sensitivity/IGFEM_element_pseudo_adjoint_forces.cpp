/*
Created by Marcus Tan on 12/23/2014
Modified on 1/6/2015
Copyright 2014 University of Illinois 
Purpose: this function calculates the IGFEM element contribution to the 
         pseudo force, adjoint force and remainder term of the objective
         function gradient
CAUTION: The extra term in DP of the pseudo force coming from 
         IGFEM has not been added for quadrilateral child element.
REMARK: i) The case of prescribed heat flux changing with design parameter has 
            been implemented here
INPUT:
    nodeCoords:
    paLocOrigNodes: local number of original nodes of the parent element wrt to itself
    elemHeatSource: element heat source
    gauss1:
    parent1:
        parent1.children(i).locPaNodes: local number of child nodes wrt parent element
        parent1.children(i).locPaEnNodes: local number of enrichment node in child wrt parent
        parent1.children(i).locEnNodes: local number of enrichment node in child wrt itself
        parent1.children(i).conductivity
	channels:
		channels.mcf: a vector of mass flow rate x heat capacity of each channel
		channels.model: only MEAN_TEMP has been implemented
    UURel: nodal solutions of the element    
OUTPUT:
    FpseudoEl: nNodesPerElem x nDesignParams matrix of the element pseudo force
               Note: the global pseudo force is of size nDesignParams x nNodes 
    FadjEl: vector of size nNodesPerElem of the element adjoint force
    gradObjRemEl: contribution of element to the remainder term in the objective function gradient
          
*/
#include "sensitivity.h"
#include "armadillo"
#include <cstddef>
#include <math.h> // pow
#include <mex.h>

#define D2MCF_TOL 1.0e-13
//#define SOURCE_FUNCTION
//#define SHOW_DK_DP

namespace igfem
{

void assemble_DK_matrix_integrand(arma::mat& DKIntegrand,
                                  double& divv,
                                  arma::vec& D2W,
                                  arma::vec& D2N,
                                  const arma::mat& vel,
                                  std::size_t nNodes,
                                  const arma::uvec& paLocOrigNodes,
                                  const arma::uvec& chLocPaEnrichNodes,
                                  const arma::vec& N,
                                  const arma::vec& W,
                                  const arma::mat& B,
                                  const arma::mat& Bel,
                                  const arma::mat& DNchEn, 
                                  const arma::mat& BchEn,
                                  const arma::mat& invJ,
                                  const arma::mat& Cmat,
                                  bool supg,
                                  const arma::mat& channelUnitVecs,
                                  const arma::mat& D2t,
                                  const arma::vec& channelMcfs,
                                  double convectCoef)
{
    D2N.zeros(N.n_elem);
    D2N(paLocOrigNodes) = Bel*vel*N(chLocPaEnrichNodes);
    DKIntegrand.zeros(nNodes,nNodes);
    arma::mat D1v = vel*BchEn;  // derivative of velocity wrt nodal position                    
    divv = arma::trace(D1v); // divergence of velocity
    arma::mat D2J = vel*DNchEn; // derivative of Jacobian wrt one design parameter  
	// derivative of the enrichment part of the B matrix wrt one design parameter	
    arma::mat D2BchEn = -DNchEn*invJ*D2J*invJ;  
    
    if (supg && channelUnitVecs.n_cols)
    {
        D2_supg_weighing_function(D2W,
                                  channelUnitVecs,
                                  D2t,
                                  chLocPaEnrichNodes,
                                  D2N,
                                  B,
                                  D2BchEn,
                                  channelMcfs);
    }
    else
        D2W = D2N;

    DKIntegrand(chLocPaEnrichNodes,chLocPaEnrichNodes) 
        = BchEn*Cmat*D2BchEn.t() + D2BchEn*Cmat*BchEn.t(); 

    DKIntegrand(paLocOrigNodes,chLocPaEnrichNodes)
        = Bel*Cmat*D2BchEn.t();

    DKIntegrand(chLocPaEnrichNodes,paLocOrigNodes)
        = D2BchEn*Cmat*Bel.t();
   
    
    DKIntegrand += convectCoef*(D2W*N.t() + W*D2N.t())
                   + (B*Cmat*B.t() + convectCoef*W*N.t())*divv;

}

void assemble_channel_DK_matrix_integrand(arma::mat& DKIntegrand,
                                          arma::mat& D2BchEn,
                                          double& surfDivv,
                                          double mcf,
                                          double D2mcf,
                                          const arma::vec& channelVec,
                                          double vecLength,
                                          const arma::mat& projMat,
                                          const arma::vec& D2t,
                                          const arma::mat& vel,
                                          std::size_t nNodes, // const arma::uvec& paLocOrigNodes,
                                          const arma::uvec& chLocPaEnrichNodes,
                                          const arma::vec& W, // weighing function
                                          const arma::mat& B, // const arma::mat& Bel,
                                          const arma::mat& DNchEn, 
                                          const arma::mat& BchEn,
                                          const arma::mat& invJ)
{
    arma::mat D1v = vel*BchEn;  // derivative of velocity wrt nodal position                    
    //double surfDivv = arma::trace(D1v*projMat); // surface divergence of velocity
    surfDivv = arma::trace(D1v*projMat); // surface divergence of velocity
    arma::mat D2J = vel*DNchEn; // derivative of Jacobian wrt one design parameter  
	// derivative of the enrichment part of the B matrix wrt one design parameter	
    D2BchEn = -DNchEn*invJ*D2J*invJ;  
    arma::vec D2mcftB = arma::zeros<arma::vec>(nNodes);

    D2mcftB(chLocPaEnrichNodes) = mcf*D2BchEn*channelVec;
    D2mcftB += B*(mcf*(D2t*vecLength + surfDivv*channelVec)+D2mcf*channelVec);

    //arma::vec D2N = arma::zeros<arma::vec>(N.n_elem);
    //D2N(paLocOrigNodes) = DNel*invJel*vel*N(chLocPaEnrichNodes);
    //D2N(paLocOrigNodes) = Bel*vel*N(chLocPaEnrichNodes);
     
    //DKIntegrand = W*D2mcftB.t() + mcf*D2W*arma::trans(B*channelVec);
    DKIntegrand = W*D2mcftB.t();
}
// DtDpts = [dtx/dx1,dtx/dy1,dtx/dx2,dtx/dy2;               
//           dty/dx1,dty/dy1,dty/dx2,dty/dy2]
void channel_unit_vec_der_wrt_end_pts(arma::mat& DtDpts,
                                      arma::uvec& channelLocPaEnNodesOffset,
                                      const arma::uvec& locPaNodes,
                                      const arma::uvec& channelLocNodes,
                                      const arma::mat& pts,
                                      double vecLength)
{
   int chanLocNode0Offset = locPaNodes(channelLocNodes(0)) - N_ORIGINAL_NODES;
   int chanLocNode1Offset = locPaNodes(channelLocNodes(1)) - N_ORIGINAL_NODES;
    
   if (chanLocNode0Offset >= 0 && chanLocNode1Offset >= 0)
   {
       DtDpts.set_size(2,4);
       DtDpts(0,0) = -pow(pts(1,0) - pts(1,1),2.0);
       DtDpts(1,0) = (pts(1,0) - pts(1,1))*(pts(0,0) - pts(0,1));
       DtDpts(0,1) = DtDpts(1,0);
       DtDpts(1,1) = -pow(pts(0,0) - pts(0,1),2.0);
       DtDpts(0,2) = -DtDpts(0,0);
       DtDpts(1,2) = -DtDpts(1,0);
       DtDpts(0,3) = -DtDpts(1,0);
       DtDpts(1,3) = -DtDpts(1,1);
       DtDpts /= pow(vecLength,3.0);
       channelLocPaEnNodesOffset.set_size(2);
       channelLocPaEnNodesOffset(0) = chanLocNode0Offset;
       channelLocPaEnNodesOffset(1) = chanLocNode1Offset;
   }
   else if (chanLocNode0Offset >= 0)
   {
       DtDpts.set_size(2,2);
       DtDpts(0,0) = -pow(pts(1,0) - pts(1,1),2.0);
       DtDpts(1,0) = (pts(1,0) - pts(1,1))*(pts(0,0) - pts(0,1));
       DtDpts(0,1) = DtDpts(1,0);
       DtDpts(1,1) = -pow(pts(0,0) - pts(0,1),2.0);
       DtDpts /= pow(vecLength,3.0);
       channelLocPaEnNodesOffset.set_size(1);
       channelLocPaEnNodesOffset(0) = chanLocNode0Offset;

   }
   else if (chanLocNode1Offset >= 0)
   {
       DtDpts.set_size(2,2);
       DtDpts(0,0) = pow(pts(1,0) - pts(1,1),2.0);
       DtDpts(1,0) = (pts(1,0) - pts(1,1))*(pts(0,1) - pts(0,0));
       DtDpts(0,1) = DtDpts(1,0);
       DtDpts(1,1) = pow(pts(0,0) - pts(0,1),2.0);
       DtDpts /= pow(vecLength,3.0);

       channelLocPaEnNodesOffset.set_size(1);
       channelLocPaEnNodesOffset(0) = chanLocNode1Offset;

   }
   else
   {
       DtDpts.set_size(0,0);
       channelLocPaEnNodesOffset.set_size(0);
   }
}

void D2_channel_unit_vec(arma::cube& D2t,
                         const arma::cube& vel,
                         const arma::uvec& locPaNodes,
                         const arma::umat& channelNodes,
                         const arma::umat& channelLocNodes,
                         const arma::mat& nodeCoords,
                         const arma::vec& vecLengths)
{
    D2t.zeros(nodeCoords.n_rows,channelNodes.n_cols,vel.n_slices);
    arma::mat DtDpts;
    arma::uvec channelLocPaEnNodesOffset;
    for (std::size_t i = 0; i < channelNodes.n_cols; i++)
    {
	    channel_unit_vec_der_wrt_end_pts(DtDpts,
									     channelLocPaEnNodesOffset,
									     locPaNodes,
									     channelLocNodes.col(i),
									     nodeCoords.cols(channelNodes.col(i)),
									     vecLengths(i));
        for (std::size_t j = 0; j < vel.n_slices; j++)
            for (std::size_t k = 0; k < channelLocPaEnNodesOffset.n_elem; k++)
            {
                D2t.slice(j).col(i) 
                    += DtDpts.cols(2*k,2*k+1)
                      *vel.slice(j).col(channelLocPaEnNodesOffset(k));
            }
    }
}

#ifdef SHOW_DK_DP
void IGFEM_element_pseudo_adjoint_forces(arma::mat& FpseudoEl,
                           arma::vec& FadjEl,
                           arma::vec& gradObjRemEl,
                           double& objValEl, 
                           arma::vec& Fadj1normEl,
                           arma::vec& grad1normRemEl,
                           double& oneNormValEl, 
                           arma::mat& DKel, arma::vec& DPel,
                           const ObjOptions& objOpt,
                           const arma::mat& nodeCoords,
                           const arma::uvec& paLocOrigNodes,
                           double elemHeatSource,
                           const Convection& convect,
                           const gauss& gauss1,
                           const parent& parent1,
                           const chanNetwork& channels,
						   const Neumann& neumann1,
                           const arma::umat& surfLocNodes,
                           const arma::vec& UURel,
                           bool supg,
                           bool calcGrad)
#else
void IGFEM_element_pseudo_adjoint_forces(arma::mat& FpseudoEl,
                           arma::vec& FadjEl,
                           arma::vec& gradObjRemEl,
                           double& objValEl, 
                           arma::vec& Fadj1normEl,
                           arma::vec& grad1normRemEl,
                           double& oneNormValEl, 
                           const ObjOptions& objOpt,
                           const arma::mat& nodeCoords,
                           const arma::uvec& paLocOrigNodes,
                           double elemHeatSource,
                           const Convection& convect,
                           const gauss& gauss1,
                           const parent& parent1,
                           const chanNetwork& channels,
						   const Neumann& neumann1,
                           const arma::umat& surfLocNodes,
                           const arma::vec& UURel,
                           bool supg,
                           bool calcGrad)

#endif
{
    #ifdef SHOW_DK_DP
        DKel.zeros(parent1.nodes.n_elem,parent1.nodes.n_elem);
        DPel.zeros(parent1.nodes.n_elem);
    #endif
    
    if (calcGrad)
    {
        FpseudoEl.zeros(parent1.nodes.n_elem,objOpt.nDesignParams);
        FadjEl.zeros(parent1.nodes.n_elem);
        gradObjRemEl.zeros(objOpt.nDesignParams);
    }
    objValEl = 0.0;

    if (objOpt.calcOneNorm) 
    {
        if (calcGrad)
        {
            Fadj1normEl.zeros(parent1.nodes.n_elem);
            grad1normRemEl.zeros(objOpt.nDesignParams);
        }
        oneNormValEl = 0.0;
    }

    // find design parameters without velocity
    // as the term D2mcf in DK may have non-zero
    // contribution even without a velocity
    arma::uvec designParamNoVelFlag = arma::ones<arma::uvec>(objOpt.nDesignParams);
    for (std::size_t j = 0; j < parent1.designParamNum.n_elem; j++)
        designParamNoVelFlag(parent1.designParamNum(j)) = 0;
    
    arma::vec N; // shape functions for IGFEM element
    arma::mat B; // B-matrix of IGFEM element
    arma::mat Bel; // B-matrix of original element
    arma::mat DNel; // DN-matrix of original element
    arma::mat Xpa = nodeCoords.cols(parent1.nodes); // all global coordinates of IGFEM element
    arma::mat Xel = Xpa.cols(paLocOrigNodes); // global coordinates of original nodes
    arma::vec Xglo(nodeCoords.n_rows); // global coordinates 
    
    arma::mat Cmat = parent1.conductivity*arma::eye<arma::mat>(
                                          nodeCoords.n_rows,
                                          nodeCoords.n_rows);

    // parent element is a linear element.
    // calcBel ensures that the B matrix is only calculated once
    bool calcBel = false; // Bel calculated here instead of in child_NDNBJ
    shape_function_2D(N,DNel,arma::zeros<arma::vec>(2));
    arma::mat Jel = Xel*DNel;
    arma::mat invJel = arma::inv(Jel);
    Bel = DNel*invJel;
	
    // intermediate variables for pseudo force
    arma::mat DKIntegrand;
    arma::vec DPIntegrand;
    arma::vec D2N = arma::zeros<arma::vec>(parent1.nodes.n_elem);
	
	// for supg
    arma::vec Bsw; // diffusive weight function for supg
    arma::vec W, D2W;

    // radiation
    bool linear = true;
    double linearRadCoef = 0.0;
    double convectRadCoef4K = 0.0; // convection coefficient + coefficient arising from radiation in the stiffness matrix 
    double threeEpsSB = 0.0;
    double threeEpsSBTsq = 0.0;
    double epsSBtimesTcube = 0.0;
 
    if (!convect.linearRad && convect.epsSB)
    {
        linear = false;
        threeEpsSB = 3.0*convect.epsSB;
        elemHeatSource += convect.coef*convect.Tamb + convect.epsSB*pow(convect.Tamb,4.0);
    }
    else
    { 
        // The linearized radiation coefficient has been added 
        // in the Mex function
        // Care must be taken to use Kelvin for the temperature
        // in the linearization
        //linearRadCoef = 4*convect.epsSB*pow(convect.Tamb,3.0);
        convectRadCoef4K = convect.coef + linearRadCoef;
        elemHeatSource += convectRadCoef4K*convect.Tamb;
    }
    
    // for body source function
    #ifdef SOURCE_FUNCTION
        double fb;
        arma::vec Dfb;
        arma::vec Nch;
    #endif
    for (std::size_t ch = 0; ch < parent1.children.n_elem; ch++)
    {
        double factor;
        double detJ;
        // if the child element is a linear element,
        // its B matrix and Jacobian are constants
        // calcBJch ensures that these quantities are calcuated
        // once only per child element
        bool calcBJch = true;
        arma::mat DNch; // DN-matrix of child element
        arma::mat Bch; // B-matrix of child element
        arma::mat invJ; // invJ is the inverse of the Jacobian of the child element
        arma::mat Xch = Xpa.cols(parent1.children(ch).locPaNodes);
        arma::vec locCoord(nodeCoords.n_rows);
        arma::uvec locPaEnNodesOffset = parent1.children(ch).locPaEnNodes 
                                       -N_ORIGINAL_NODES;
        arma::vec v; // velocity at a point
        arma::mat D1v; // derivative of velocity wrt global coordinates
        double divv; // divergence of velocity
        double uh; // finite element temperature
        double uhloffset; // finite element temperature - offset
        // for channel contribution
		arma::mat channelVecs, channelUnitVecs;
       	arma::vec vecLengths;        
        arma::cube D2t;
        arma::mat D2BchEn;

        if (parent1.children(ch).channelNum.n_elem > 0)
        {      
            channelVecs.set_size(2,parent1.children(ch).channelNum.n_elem);
            channelUnitVecs.set_size(channelVecs.n_rows,channelVecs.n_cols);
            vecLengths.set_size(channelVecs.n_cols);
            for (std::size_t chanNum = 0; chanNum < parent1.children(ch).channelNum.n_elem; chanNum++)
            {
                channelVecs.col(chanNum)
                        = nodeCoords.col(parent1.children(ch).channelNodes(1,chanNum))
                         -nodeCoords.col(parent1.children(ch).channelNodes(0,chanNum));
                vecLengths(chanNum) = arma::norm(channelVecs.col(chanNum));
                channelUnitVecs.col(chanNum) = channelVecs.col(chanNum)
                                                /vecLengths(chanNum);
            }
            
            D2_channel_unit_vec(D2t,
                                parent1.vel,
                                parent1.children(ch).locPaNodes,
                                parent1.children(ch).channelNodes,
                                parent1.children(ch).channelLocNodes,
                                nodeCoords,
                                vecLengths);
        
        }
        else
            D2t.zeros(1,1,parent1.vel.n_slices);
                   
        if (parent1.children(ch).shape == TRIANGLE)
        {
            for (std::size_t i = 0; i < gauss1.elem.n_cols; i++)
            {               
                child_element_NDNBJ(  N,
                                      B,
                                      DNch,
                                      Bch,
                                      DNel,
                                      Bel,
                                      invJ,
                                      detJ,
                                      calcBJch,
                                      calcBel,
                                      Xglo,
                                      parent1.nodes.n_elem,
                                      paLocOrigNodes,
                                      parent1.children(ch).locPaEnNodes,
                                      parent1.children(ch).locEnNodes,
                                      Xel,
                                      Xch,
                                      gauss1.elem(arma::span(0,gauss1.elem.n_rows-2),i),
                                      parent1.children(ch).shape);
				
                if (supg  && parent1.children(ch).channelNum.n_elem)
                    supg_weighing_function(W,channelUnitVecs,N,B,
                                           channels.mcf(parent1.children(ch).channelNum));
                else
                    W = N;
			    
                factor = detJ*gauss1.elem(gauss1.elem.n_rows-1,i);
               
                uh = arma::dot(N,UURel);
                if (!linear)
                {
                    epsSBtimesTcube = convect.epsSB*pow(uh,3.0);
                    threeEpsSBTsq = threeEpsSB*pow(uh,2.0);
                }
                uhloffset = uh - objOpt.offset;
                if (objOpt.intDomain == WHOLE)
                {
                    objValEl += pow(uhloffset,objOpt.normp)*factor;
                    if (calcGrad) FadjEl += N*objOpt.normp*pow(uhloffset,objOpt.normp-1)*factor;
                }
                if (objOpt.calcOneNorm)
                {
                    oneNormValEl += uh*factor;
                    if (calcGrad) Fadj1normEl += N*factor;
                }
               
                #ifdef SOURCE_FUNCTION                 
                    body_source(fb,Dfb,Xglo);
                #endif
                
                if (calcGrad)
                {
                    if (!linear)
                        convectRadCoef4K = convect.coef + epsSBtimesTcube;

                    for (std::size_t j = 0; j < parent1.vel.n_slices; j++)
                    {
                        assemble_DK_matrix_integrand(DKIntegrand,
                                                     divv,
                                                     D2W,
                                                     D2N,
                                                     parent1.vel.slice(j).cols(locPaEnNodesOffset),
                                                     parent1.nodes.n_elem,
                                                     paLocOrigNodes,
                                                     parent1.children(ch).locPaEnNodes,
                                                     N,
                                                     W,
                                                     B,
                                                     Bel,
                                                     DNch.rows(parent1.children(ch).locEnNodes), 
                                                     Bch.rows(parent1.children(ch).locEnNodes),
                                                     invJ,
                                                     Cmat,
                                                     supg,
                                                     channelUnitVecs,
                                                     D2t.slice(j),
                                                     channels.mcf(parent1.children(ch).channelNum),
                                                     convectRadCoef4K);
                        if (!linear)
                            DKIntegrand += threeEpsSBTsq*arma::dot(UURel,D2N)*W*N.t();
                        
                        #ifdef SOURCE_FUNCTION       
                            v = parent1.vel.slice(j).cols(locPaEnNodesOffset)*N(parent1.children(ch).locPaEnNodes);                       
                            DPIntegrand = (W*divv+D2W)*(elemHeatSource + fb) + W*Dfb.t()*v;
                        #else
                            DPIntegrand = (W*divv+D2W)*elemHeatSource;                    
                        #endif
                        //DPIntegrand -= B*Cmat*arma::trans(B*D1v)*UPel; // should be wrong
                        #ifdef SHOW_DK_DP       
                            DKel += DKIntegrand*factor; 
                            DPel += DPIntegrand*factor;
                        #endif
                     
                        FpseudoEl.col(parent1.designParamNum(j)) 
                                += (DPIntegrand - DKIntegrand*UURel)*factor;
                        if (objOpt.calcOneNorm) 
                        {
                            grad1normRemEl(parent1.designParamNum(j))                         
                                += (uh*divv + arma::dot(UURel,D2N))*factor;                            
                        }       

                        if (objOpt.intDomain == WHOLE)
                        {
                            gradObjRemEl(parent1.designParamNum(j))
                                += (pow(uhloffset,objOpt.normp)*divv+objOpt.normp*pow(uhloffset,objOpt.normp-1)*arma::dot(UURel,D2N))*factor;                            
                            //gradObjRemEl(parent1.designParamNum(j)) 
                            //    += objOpt.normp*pow(uh2,objOpt.normp-1)*arma::dot(UPel,B*v)*factor;
                        }
                    }
                }
           
            }
        }
        else if (parent1.children(ch).shape == QUADRILATERAL)
        {
          
            for (std::size_t i = 0; i < gauss1.quadElem.n_cols; i++)
            {
         
                child_element_NDNBJ(  N,
                                      B,
                                      DNch,
                                      Bch,
                                      DNel,
                                      Bel,
                                      invJ,
                                      detJ,
                                      calcBJch,
                                      calcBel,
                                      Xglo,
                                      parent1.nodes.n_elem,
                                      paLocOrigNodes,
                                      parent1.children(ch).locPaEnNodes,
                                      parent1.children(ch).locEnNodes,
                                      Xel,
                                      Xch,
                                      gauss1.elem(arma::span(0,gauss1.elem.n_rows-2),i),
                                      parent1.children(ch).shape);
                
                if (supg  && parent1.children(ch).channelNum.n_elem)
                    supg_weighing_function(W,channelUnitVecs,N,B,
                                           channels.mcf(parent1.children(ch).channelNum));
                else
                    W = N;

                factor = detJ*gauss1.quadElem(gauss1.quadElem.n_rows-1,i);
                
                uh = arma::dot(N,UURel) - objOpt.offset;
                if (!linear)
                {
                    epsSBtimesTcube = convect.epsSB*pow(uh,3.0);
                    threeEpsSBTsq = threeEpsSB*pow(uh,2.0);
                }
                uhloffset = uh - objOpt.offset;
                if (objOpt.intDomain == WHOLE)
                {
                    objValEl += pow(uhloffset,objOpt.normp)*factor;
                    if (calcGrad) FadjEl += N*objOpt.normp*pow(uhloffset,objOpt.normp-1)*factor;
                }

                if (objOpt.calcOneNorm) 
                {
                    oneNormValEl += uh*factor;
                    if (calcGrad) Fadj1normEl += N*factor;
                }
                #ifdef SOURCE_FUNCTION
                    body_source(fb,Dfb,Xglo);
                #endif
                
                if (calcGrad)
                {
                    if (!linear)
                        convectRadCoef4K = convect.coef + epsSBtimesTcube;
                       
                    for (std::size_t j = 0; j < parent1.vel.n_slices; j++)
                    {
                        assemble_DK_matrix_integrand(DKIntegrand,
                                                     divv,
                                                     D2W,
                                                     D2N,
                                                     parent1.vel.slice(j).cols(locPaEnNodesOffset),
                                                     parent1.nodes.n_elem,
                                                     paLocOrigNodes,
                                                     parent1.children(ch).locPaEnNodes,
                                                     N,
                                                     W,
                                                     B,
                                                     Bel,
                                                     DNch.rows(parent1.children(ch).locEnNodes), 
                                                     Bch.rows(parent1.children(ch).locEnNodes),
                                                     invJ,
                                                     Cmat,
                                                     supg,
                                                     channelUnitVecs,
                                                     D2t.slice(j),
                                                     channels.mcf(parent1.children(ch).channelNum),
                                                     convectRadCoef4K);
                        if (!linear)
                            DKIntegrand += threeEpsSBTsq*arma::dot(UURel,D2N)*W*N.t();
                        
                        #ifdef SOURCE_FUNCTION
                            v = parent1.vel.slice(j).cols(locPaEnNodesOffset)*N(parent1.children(ch).locPaEnNodes);       
                            DPIntegrand = (W*divv+D2W)*(elemHeatSource + fb)
                                          +W*Dfb.t()*v;
                        #else
                            DPIntegrand = (W*divv+D2W)*elemHeatSource;
                        #endif
                                            
                        #ifdef SHOW_DK_DP       
                            DKel += DKIntegrand*factor; 
                            DPel += DPIntegrand*factor;
                        #endif
                        
                        FpseudoEl.col(parent1.designParamNum(j)) 
                                += (DPIntegrand - DKIntegrand*UURel)*factor;          
                        if (objOpt.calcOneNorm)
                        {
                            grad1normRemEl(parent1.designParamNum(j))                         
                                += (uh*divv + arma::dot(UURel,D2N))*factor;                            
                        }
                        if (objOpt.intDomain == WHOLE)
                        {
                            gradObjRemEl(parent1.designParamNum(j))
                                += (pow(uhloffset,objOpt.normp)*divv+objOpt.normp*pow(uhloffset,objOpt.normp-1)*arma::dot(UURel,D2N))*factor;
                        }
                    }
                }
            }
        }
        else
        {
            //std::cerr << "IGFEM_element_pseudo_adjoint_forces:" 
            //          << "unknown child element shape" << std::endl;
            mexErrMsgIdAndTxt("IGFEM_element_pseudo_adjoint_forces:child_shape",
                              "IGFEM_element_pseudo_adjoint_forces:unknown child element shape\n");
        }
	    ////////////////////////////////////////////////////////////////
        // NOTE: The loop is continued to the next child element      //
        //       if the sensitivity is not needed and the objective   //
        //       function is not integrated over the channels!        // 
        ////////////////////////////////////////////////////////////////
        if (!calcGrad && objOpt.intDomain != CHANNEL) continue;
  
		if (channels.model == MEAN_TEMP)
        {
            arma::vec unitNormal(nodeCoords.n_rows);
            arma::uvec channelLocPaEnNodesOffset;
            arma::mat projMat, DtDpts;
            int edge;
            double surfDivv;
            double p5gaussWeight;
            factor = 0.0;
            uh = 0.0;
            for (std::size_t chanNum = 0; chanNum < parent1.children(ch).channelNum.n_elem; chanNum++)
            {
                if(fabs(channels.mcf(parent1.children(ch).channelNum(chanNum))) < MCFTOL)
                    continue;
                unitNormal(0) = -channelUnitVecs(1,chanNum);
                unitNormal(1) =  channelUnitVecs(0,chanNum);
                projMat = arma::eye<arma::mat>(nodeCoords.n_rows,nodeCoords.n_rows)
                         - unitNormal*unitNormal.t();    

                edge = edge_number(parent1.children(ch).channelLocNodes.col(chanNum),
                                   parent1.children(ch).shape);    
                
				for (std::size_t i = 0; i < gauss1.line.n_cols; i++)
				{
					Xglo = nodeCoords.col(parent1.children(ch).channelNodes(0,chanNum)) 
						  + gauss1.line(0,i)*channelVecs.col(chanNum);
					locCoord = local_coord_2D_along_edge(Xglo,
														 Xch, 
				  									     parent1.children(ch).shape,
														 edge);
					
                    child_element_NDNBJ(  N,
										  B,
										  DNch,
										  Bch,
										  DNel,
										  Bel,
										  invJ,
										  detJ,
										  calcBJch,
										  calcBel,
                                          Xglo,
										  parent1.nodes.n_elem,
										  paLocOrigNodes,
										  parent1.children(ch).locPaEnNodes,
										  parent1.children(ch).locEnNodes,
										  Xel,
										  Xch,
										  locCoord,
										  parent1.children(ch).shape);
                    p5gaussWeight = 0.5*gauss1.line(1,i);
                    if (objOpt.intDomain == CHANNEL)
                    {
                        factor = vecLengths(chanNum)*p5gaussWeight;
                        uhloffset = arma::dot(N,UURel) - objOpt.offset;
                        // I assume that the channels are never at the boundary of the domain
                        // So they are always shared by two child elements. Hence 0.5
                        objValEl += pow(uhloffset,objOpt.normp)*factor;
                        //objValEl += factor; // length for debugging
                        if (calcGrad) FadjEl += N*objOpt.normp*pow(uhloffset,objOpt.normp-1)*factor;
                    }
                    if (calcGrad)
                    {
                        if (supg  && parent1.children(ch).channelNum.n_elem)
                            supg_weighing_function(W,channelUnitVecs,N,B,
                                                   channels.mcf(parent1.children(ch).channelNum));
                        else
                            W = N;
                            
                        for (std::size_t j = 0; j < parent1.vel.n_slices; j++)
                        {
                            assemble_channel_DK_matrix_integrand(DKIntegrand,
                                                                 D2BchEn,
                                                                 surfDivv,
                                                                 channels.mcf(parent1.children(ch).channelNum(chanNum)),
                                                                 channels.D2mcf(parent1.designParamNum(j),
                                                                                parent1.children(ch).channelNum(chanNum)),
                                                                 channelVecs.col(chanNum),
                                                                 vecLengths(chanNum),
                                                                 projMat,
                                                                 D2t.slice(j).col(chanNum),
                                                                 parent1.vel.slice(j).cols(locPaEnNodesOffset),
                                                                 parent1.nodes.n_elem, //paLocOrigNodes,
                                                                 parent1.children(ch).locPaEnNodes,
                                                                 W,
                                                                 B, // Bel,
                                                                 DNch.rows(parent1.children(ch).locEnNodes), 
                                                                 Bch.rows(parent1.children(ch).locEnNodes),
                                                                 invJ);

                            v = parent1.vel.slice(j).cols(locPaEnNodesOffset)*N(parent1.children(ch).locPaEnNodes);                   
                            D2N(paLocOrigNodes) = Bel*v;
                            
                            if (supg && parent1.children(ch).channelNum.n_elem)
                                D2_supg_weighing_function(D2W,
                                                          channelUnitVecs,
                                                          D2t.slice(j),
                                                          parent1.children(ch).locPaEnNodes,
                                                          D2N,
                                                          B,
                                                          D2BchEn,
                                                          channels.mcf(parent1.children(ch).channelNum));

                            else
                                D2W = D2N;
                            
                            DKIntegrand += channels.mcf(parent1.children(ch).channelNum(chanNum))
                                         *D2W*arma::trans(B*channelVecs.col(chanNum));

                            #ifdef SHOW_DK_DP             
                                DKel += DKIntegrand*p5gaussWeight;
                            #endif
                            //DPIntegrand = -N*arma::trans(B*D1v*channels.mcf(parent1.children(ch).channelNum(chanNum))
                            //                *channelVec)*UPel; // should be wrong
                            FpseudoEl.col(parent1.designParamNum(j)) -= DKIntegrand*UURel*p5gaussWeight;
                      
                            if (objOpt.intDomain == CHANNEL)
                            {
                                // Ensure that factor = 0.5*vecLengths(chanNum)*gauss1.line(1,i); 
                                factor = vecLengths(chanNum)*0.5*gauss1.line(1,i);
                                gradObjRemEl(parent1.designParamNum(j))
                                    += (pow(uhloffset,objOpt.normp)*surfDivv+objOpt.normp*pow(uhloffset,objOpt.normp-1)*arma::dot(UURel,D2N))*factor;
                            }
                        }
                    }
                    // reuse N for another purpose here. don't confuse with the shape function.
                    N = -0.5*W*arma::trans(B*channelVecs.col(chanNum))*UURel*gauss1.line(1,i);
                    for (size_t j = 0; j < objOpt.nDesignParams; j++)
                    {
                        
                        if (designParamNoVelFlag(j) 
                            && fabs(channels.D2mcf(j,parent1.children(ch).channelNum(chanNum))) > D2MCF_TOL)
                        {
                            FpseudoEl.col(j)
                                += channels.D2mcf(j,parent1.children(ch).channelNum(chanNum))*N;
                        }
                    }					                        
				}
            }
        }
		else if (channels.model == CONST_HEAT)
		{
			//std::cerr << "IGFEM_element_pseudo_adjoint_forces:" 
            //          << "const heat channel model not yet implemented" << std::endl;
            mexErrMsgIdAndTxt("IGFEM_element_pseudo_adjoint_forces:const_heat",
                              "IGFEM_element_pseudo_adjoint_forces:const heat channel model not yet implemented\n");
		}
		else
		{
            //std::cerr << "IGFEM_element_pseudo_adjoint_forces:" 
            //          << "unknown channel model" << std::endl;	
            mexErrMsgIdAndTxt("IGFEM_element_pseudo_adjoint_forces:unknown_model",
                              "IGFEM_element_pseudo_adjoint_forces:unknown model\n");

		}

        // Assumption: there is no enrichment node at the boundary of the domain
        if (neumann1.surf >= 0 && neumann1.surf < 3)
        {
            if (arma::any(parent1.children(ch).locPaNodes == surfLocNodes(0,neumann1.surf)) 
                && arma::any(parent1.children(ch).locPaNodes == surfLocNodes(1,neumann1.surf)))
            {
                double dS = 0;
                if (neumann1.surf == 0)
                    dS = sqrt(Jel(0,0)*Jel(0,0) + Jel(1,0)*Jel(1,0));
                else if (neumann1.surf == 1)
                    dS = sqrt(pow(Jel(0,1) - Jel(0,0),2) + pow(Jel(1,1) - Jel(1,0),2));
                else if (neumann1.surf == 2)
                    dS = sqrt(Jel(0,1)*Jel(0,1) + Jel(1,1)*Jel(1,1));

				for (std::size_t i = 0; i < gauss1.line.n_cols; i++)
				{
		            if (neumann1.surf == 0)
                    {
                        locCoord(0) = gauss1.line(0,i);
                        locCoord(1) = 0;
                    }
                    else if (neumann1.surf == 1)
                    {
                        locCoord(0) = gauss1.line(0,i);
                        locCoord(1) = 1.0 - gauss1.line(0,i);
                    }
                    else if (neumann1.surf == 2)
                    {
                        locCoord(0) = 0;
                        locCoord(1) = gauss1.line(0,i);
                    }
                    child_element_NDNBJ(  N,
										  B,
										  DNch,
										  Bch,
										  DNel,
										  Bel,
										  invJ,
										  detJ,
										  calcBJch,
										  calcBel,
                                          Xglo,
										  parent1.nodes.n_elem,
										  paLocOrigNodes,
										  parent1.children(ch).locPaEnNodes,
										  parent1.children(ch).locEnNodes,
										  Xel,
										  Xch,
										  locCoord,
										  parent1.children(ch).shape);

					for (std::size_t j = 0; j < parent1.vel.n_slices; j++)
					{
                        v = parent1.vel.slice(j).cols(locPaEnNodesOffset)*N(parent1.children(ch).locPaEnNodes);                   
                        D2N(paLocOrigNodes) = Bel*v;

                        // Note: heat flux term does not require modified weight function
                        /*
                        if (supg && parent1.children(ch).channelNum.n_elem)
                            D2_supg_weighing_function(D2W,
                                                      channelUnitVecs,
                                                      D2t.slice(j),
                                                      parent1.children(ch).locPaEnNodes,
                                                      D2N,
                                                      B,
                                                      D2BchEn,
                                                      channels.mcf(parent1.children(ch).channelNum));

                        else
                            D2W = D2N;
                        */
                        divv = arma::trace(parent1.vel.slice(j).cols(locPaEnNodesOffset)
                                           *Bch.rows(parent1.children(ch).locEnNodes));
						// FpseudoEl.col(parent1.designParamNum(j)) 
                        //    += (N*divv + D2W)*dS*gauss1.line(1,i); 
						FpseudoEl.col(parent1.designParamNum(j)) 
                            += (N*divv + D2N)*dS*gauss1.line(1,i); 

                    }


                }
                
            }
        }
    } 
    
    /*
    for (std::size_t j = 0; j < parent1.vel.n_slices; j++)
    {
            FpseudoEl.col(parent1.designParamNum(j)) 
                = - DK.slice(j)*UURel + DP.col(j);
    }
    */
}   


}

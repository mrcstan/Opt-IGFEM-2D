/*
Created by Marcus Tan on 12/23/2014
Modified on 2/21/2015
Copyright 2014 University of Illinois 
Purpose: this function assembles the pseudo forces for the direct method,
         the RHS for the adjoint equation and an additional term for computation
         of the derivative of the objective function
Formulation:
    Direct method:  K*U_star = trans(Fpseudo) 
                    gradObj = trans(U_star)*Fadj + gradObjRem
    Adjoint method: K*Lambda = Fadj
                    gradObj = Fpseudo*Lambda + gradObjRem

INPUT: 

OUTPUT:
    Fpseudo: a nDesignParams x nNodes pseudo force matrix 
             (also called trans(Pps) in notes)
             for column-wise write access, I defined it 
             that way instead of the transpose
    Fadj: an adjoint force vector of size nNodes 
          (also called trans(D1_PI) in notes)
    gradObjRem: a vector of size nDesignParams of the 
                remainder term of the gradient of the 
                objective function wrt design parameters 
                (also called D3_PI in notes)

*/

#include "sensitivity.h"
#include <omp.h>
#include <cstddef> // NULL, std::size_t
#include <iostream>
#include "armadillo"
#include "mex.h"

//#define SHOW_DK_DP

namespace igfem
{
void assemble_pseudo_adjoint_forces(arma::mat& Fpseudo,
                                    arma::vec& Fadj,
                                    arma::vec& gradObjRem,
                                    double& objVal,
                                    arma::vec& Fadj1norm,
                                    arma::vec& grad1normRem,
                                    double& oneNormVal,
                                    const ObjOptions& objOpt,
                                    const arma::mat& nodeCoords,
                                    const arma::umat& elemNodes,
                                    const arma::vec& elemHeatSource,
                                    const Convection& convect,
                                    const gauss& gauss1,
                                    const arma::field<parent>& parents,
                                    const chanNetwork& channels,
                                    const arma::field<Neumann>& elemNeumann,
                                    const arma::vec& UUR,
                                    bool supg,
                                    bool calcGrad) 
{
    // initialization
    if (calcGrad)
    {
        Fpseudo.zeros(objOpt.nDesignParams,UUR.n_elem);    
        Fadj.zeros(UUR.n_elem);
        gradObjRem.zeros(objOpt.nDesignParams);
    }
    objVal = 0.0;    
        
    if (objOpt.calcOneNorm)
    {
        oneNormVal = 0.0; 
        if (calcGrad)
        {
            Fadj1norm.zeros(UUR.n_elem);
            grad1normRem.zeros(objOpt.nDesignParams);
        }
    }
    
    
    #ifdef SHOW_DK_DP
        arma::mat DK;
        arma::vec DP;
        if (calcGrad)
        {   
            DK.zeros(UUR.n_elem,UUR.n_elem);
            DP.zeros(UUR.n_elem);
        }
        //std::cout.precision(12);
        //std::cout << std::scientific;
    #endif

    // assume that the parent element nodes are arranged such that the first nodes
    // are the original nodes
    arma::uvec paLocOrigNodes(N_ORIGINAL_NODES);
    for (std::size_t i = 0; i < N_ORIGINAL_NODES; i++)
        paLocOrigNodes(i) = i;
   
    // initialize the local node numbers corresponding to a Neumann surface
    arma::umat surfLocNodes(2,3);
    surfLocNodes << 0 << 1 << 2 << arma::endr
                 << 1 << 2 << 0 << arma::endr;

    //omp_lock_t writelock;
    //omp_init_lock(&writelock);
    omp_set_num_threads(1); // comment this line if multiple thread is required 
    #pragma omp parallel  
    {

        int threadNum = omp_get_thread_num();
        int numThreads = omp_get_num_threads();
        //int errFlag;
        
        arma::mat FpseudoEl;
        arma::vec FadjEl;
        arma::vec gradObjRemEl; 
        double objValEl;
        arma::vec Fadj1normEl;
        arma::vec grad1normRemEl; 
        double oneNormValEl;
        
        arma::mat DKel;
        arma::vec DPel;

        // calculation of element stiffness matrix 
        // and load vector is distributed to multiple threads
        int istart = threadNum*parents.n_elem/numThreads;
        int iend;
        if (threadNum == numThreads-1)
            iend = (int) parents.n_elem;
        else
            iend = istart + parents.n_elem/numThreads;
       
        // should just loop through list of IGFEM elements to prevent 
        // load inbalance for multi-threads. but right now, 
        // this doesn't matter much since code is mainly run sequentially.
        for(int i = istart; i < iend; ++i)   
        {
            //mexPrintf("Regular element %i\n",i);
            if (parents(i).type == REGULAR)
            {
                regular_element_adjoint_forces(FadjEl,
                                              objValEl,
                                              Fadj1normEl,
                                              oneNormValEl,
                                              objOpt,
                                              nodeCoords.cols(elemNodes.col(i)),
                                              elemHeatSource(i),
                                              gauss1,
                                              parents(i),
                                              channels,
                                              UUR(elemNodes.col(i)),
                                              calcGrad);
                
                /* uncomment this when the pseudo force is non-zero 
                #pragma omp critical (pseudo)
                {
                    // the element pseudo-force matrix is of size
                    // nNodesPerElem x nDesignParams
                    Fpseudo.cols(elemNodes.col(i)) += FpseudoEl.t();
                }
                */
                if (calcGrad) 
                {
                    #pragma omp critical (adjoint)
                    {
                        Fadj(elemNodes.col(i)) += FadjEl; 
                    }
                    if (objOpt.calcOneNorm) 
                    {
                        #pragma omp critical (adjoint1norm)
                        {
                            Fadj1norm(elemNodes.col(i)) += Fadj1normEl; 
                        }
                    }
                }
                /* uncomment this when D3PI is non-zero
                #pragma omp critical (gradObjRem)
                {
                    gradObjRem += gradObjRemEl;
                }
                */
                 #pragma omp atomic
                    objVal += objValEl;
                
                if (objOpt.calcOneNorm) {
                    #pragma omp atomic
                        oneNormVal += oneNormValEl;
                }
            }
            else if (parents(i).type == IGFEM)
            {   
                //mexPrintf("IGFEM element %i\n",i);
                #ifdef SHOW_DK_DP
                    IGFEM_element_pseudo_adjoint_forces(FpseudoEl,
                                                    FadjEl,
                                                    gradObjRemEl,
                                                    objValEl, 
                                                    Fadj1normEl,
                                                    grad1normRemEl,
                                                    oneNormValEl,
                                                    DKel, DPel,
                                                    objOpt,
                                                    nodeCoords,
                                                    paLocOrigNodes,
                                                    elemHeatSource(i),
                                                    convect,
                                                    gauss1,
                                                    parents(i),
                                                    channels,
                                                    elemNeumann(i),
                                                    surfLocNodes,
                                                    UUR(parents(i).nodes),
                                                    supg,
                                                    calcGrad);
               #else
                    IGFEM_element_pseudo_adjoint_forces(FpseudoEl,
                                                    FadjEl,
                                                    gradObjRemEl,
                                                    objValEl, 
                                                    Fadj1normEl,
                                                    grad1normRemEl,
                                                    oneNormValEl,
                                                    objOpt,
                                                    nodeCoords,
                                                    paLocOrigNodes,
                                                    elemHeatSource(i),
                                                    convect,
                                                    gauss1,
                                                    parents(i),
                                                    channels,
                                                    elemNeumann(i),
                                                    surfLocNodes,
                                                    UUR(parents(i).nodes),
                                                    supg,
                                                    calcGrad);

               #endif
                                                    //UP(parents(i).nodes)
                // better to assemble the matrix and vectors for each thread
                // and invoke critical after each thread finish its computation
                // to add the thread matrix and vectors to the shared ones
                // I am not implementing this because this code is mainly
                // run sequentially.
                if (calcGrad)
                {
                    #pragma omp critical (pseudo)
                    {
                        // the element pseudo-force matrix is of size
                        // nNodesPerElem x nDesignParams
                        Fpseudo.cols(parents(i).nodes) += FpseudoEl.t();
                        //Fpseudo.rows(parents(i).nodes) += FpseudoEl;
                    }
                    #pragma omp critical (adjoint)
                    {
                        Fadj(parents(i).nodes) += FadjEl; 
                    }                
                    #pragma omp critical (gradObjRem)
                    {
                        gradObjRem += gradObjRemEl;
                    }
                    #ifdef SHOW_DK_DP
                        #pragma omp critical (DK)
                        {
                            DK(parents(i).nodes,parents(i).nodes) += DKel;
                            DP(parents(i).nodes) += DPel;
                        }
                    #endif
                }
                if (objOpt.calcOneNorm && calcGrad) {
                    #pragma omp critical (grad1normRem) 
                    {                        
                        grad1normRem += grad1normRemEl;
                    }
                    #pragma omp critical (adjoint1norm)
                    {
                        Fadj1norm(parents(i).nodes) += Fadj1normEl; 
                    }               
                }
                
                #pragma omp atomic
                    objVal += objValEl;
                
                if (objOpt.calcOneNorm) {
                    #pragma omp atomic
                        oneNormVal += oneNormValEl;
                }
               
            }
            else
                //std::cerr << "assemble_pseudo_adjoint_forces: unknown element type" 
                //          << std::endl;
                mexErrMsgIdAndTxt("assemble_pseudo_adjoint_forces:unknown_elem_type",
                                  "assemble_pseudo_adjoint_forces unknown element type");
        }        
      
    }

    #ifdef SHOW_DK_DP
        //DK.raw_print(std::cout,"DK = ");
        //DP.raw_print(std::cout,"DP = ");
        arma_mexprintf("DK = ",DK);
        arma_mexprintf("DP = ",DP);
    #endif

}

}

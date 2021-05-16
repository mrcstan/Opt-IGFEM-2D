/*
Created by Marcus Tan on 12/23/2014
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

//#define SHOW_DK_DP

namespace igfem
{
void assemble_pseudo_adjoint_forces(arma::mat& Fpseudo,
                                    arma::vec& Fadj,
                                    arma::vec& gradObjRem,
                                    double& objVal,
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
                                    bool supg) 
                                    // const arma::vec& UP
{
    // initialization
    Fpseudo.zeros(objOpt.nDesignParams,nodeCoords.n_cols);    
    Fadj.zeros(nodeCoords.n_cols);
    gradObjRem.zeros(objOpt.nDesignParams);
    objVal = 0.0;    
    
    #ifdef SHOW_DK_DP
        arma::mat DK = arma::zeros<arma::mat>(nodeCoords.n_cols,nodeCoords.n_cols);
        arma::vec DP = arma::zeros<arma::vec>(nodeCoords.n_cols);
        std::cout.precision(12);
        std::cout << std::scientific;
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
    #pragma omp parallel shared(nodeCoords,paLocOrigNodes,gauss1) 
    {

        int threadNum = omp_get_thread_num();
        int numThreads = omp_get_num_threads();
        //int errFlag;
        
        arma::mat FpseudoEl;
        arma::vec FadjEl;
        arma::vec gradObjRemEl; 
        double objValEl;
        
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
            if (parents(i).type == REGULAR)
            {
                regular_element_pseudo_adjoint_forces(FpseudoEl,
                                                      FadjEl,
                                                      gradObjRemEl,
                                                      objValEl,
                                                      objOpt,
                                                      nodeCoords.cols(elemNodes.col(i)),
                                                      elemHeatSource(i),
                                                      gauss1,
                                                      UUR(elemNodes.col(i)));
                
                /* uncomment this when the pseudo force is non-zero 
                #pragma omp critical (pseudo)
                {
                    // the element pseudo-force matrix is of size
                    // nNodesPerElem x nDesignParams
                    Fpseudo.cols(elemNodes.col(i)) += FpseudoEl.t();
                }
                */
                #pragma omp critical (adjoint)
                {
                    Fadj(elemNodes.col(i)) += FadjEl; 
                }                
                /* uncomment this when D3PI is non-zero
                #pragma omp critical (gradObjRem)
                {
                    gradObjRem += gradObjRemEl;
                }
                */
            }
            else if (parents(i).type == IGFEM)
            {                
                #ifdef SHOW_DK_DP
                    IGFEM_element_pseudo_adjoint_forces(FpseudoEl,
                                                    FadjEl,
                                                    gradObjRemEl,
                                                    objValEl, 
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
                                                    supg);
               #else
                    IGFEM_element_pseudo_adjoint_forces(FpseudoEl,
                                                    FadjEl,
                                                    gradObjRemEl,
                                                    objValEl, 
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
                                                    supg);

               #endif
                                                    //UP(parents(i).nodes)
                // better to assemble the matrix and vectors for each thread
                // and invoke critical after each thread finish its computation
                // to add the thread matrix and vectors to the shared ones
                // I am not implementing this because this code is mainly
                // run sequentially.
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
            else
                std::cerr << "assemble_pseudo_adjoint_forces: unknown element type" 
                          << std::endl;

            #pragma omp atomic
                objVal += objValEl;

        }        
      
    }

    #ifdef SHOW_DK_DP
        DK.raw_print(std::cout,"DK = ");
        DP.raw_print(std::cout,"DP = ");
    #endif

}

}

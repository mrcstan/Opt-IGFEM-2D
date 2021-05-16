/*
Created by Marcus Tan on 12/23/2014
Copyright 2014 University of Illinois 
Purpose: this function calculates the regular element contribution to the 
         pseudo force, adjoint force and remainder term of the objective
         function gradient
REMARK: i) if the conductivity does not explicityly depend on design parameters, 
        the pseudo force is identically zero.
        if the objective function and conductivity do not explicitly 
        depend on design parameters, the gradObjRemEl (D3PI) term is
        identically zero.
        however, the adjoint force term is not zero in general
        ii) The case of changing prescribed heat flux with design parameters
            has not been implemented here.
            In this case, the original instead of the modified weight function
            (original+SUPG) should be used 
INPUT:
    Xel:
    objOpt: struct variable containing info about the objective function type 
            and other relevant info (see sensitivty.h)
    gauss1:
    UURel: nodal solutions of the element   
    conductivity: (don't need this right now)

OUTPUT:
    FpseudoEl: nNodesPerElem x nDesignParams matrix of the element pseudo force
               Note: the global pseudo force is of size nDesignParams x nNodes 
    FadjEl: vector of size nNodesPerElem of the element adjoint force
    gradObjRemEl: contribution of element to the remainder term in the objective function gradient
          
*/
#include "sensitivity.h"
#include "armadillo"
#include <cstddef>
#include <iostream>
#include <math.h> // pow

namespace igfem
{

void regular_element_adjoint_forces(arma::vec& FadjEl,
                                   double& objValEl,
                                   arma::vec& Fadj1normEl,
                                   double& oneNormValEl,
                                   const ObjOptions& objOpt,
                                   const arma::mat& Xel,
                                   double elemHeatSource,
                                   const gauss& gauss1,
                                   const parent& parent1,
                                   const chanNetwork& channels,
                                   const arma::vec& UURel,
                                   bool calcGrad)
{
    if (calcGrad) FadjEl.zeros(Xel.n_cols);
    objValEl = 0.0;

    if (objOpt.calcOneNorm) 
    {
        if (calcGrad) Fadj1normEl.zeros(Xel.n_cols);
        oneNormValEl = 0.0;
    }
    arma::vec N; // shape functions for IGFEM element
    arma::mat DN; // B-matrix of IGFEM element
    arma::vec locCoord = arma::zeros<arma::vec>(Xel.n_rows); // local coordinates 
    
    shape_function_2D(N,DN,locCoord);

    arma::mat J = Xel*DN;
    double detJ = arma::det(J);
    double factor, uh, uhloffset;
    
    if (objOpt.intDomain == WHOLE)
    {
        for (std::size_t i = 0; i < gauss1.elem.n_cols; i++)
        {      
            shape_function_2D(N,DN,gauss1.elem(arma::span(0,gauss1.elem.n_rows-2),i));
            factor = detJ*gauss1.elem(gauss1.elem.n_rows-1,i);
            
            uh = arma::dot(N,UURel); 
            
            if (objOpt.calcOneNorm) 
            {
                oneNormValEl += uh*factor;
                if (calcGrad) Fadj1normEl += N*factor;
            }
            
            uhloffset = uh - objOpt.offset;
            objValEl += pow(uhloffset,objOpt.normp)*factor;
            if (calcGrad) FadjEl += N*objOpt.normp*pow(uhloffset,objOpt.normp-1)*factor;
        }
    }
    else if (objOpt.intDomain == CHANNEL)
    {
        // Note that when a channel coincides with an edge of a regular element, 
        // the derivatives are undefined. So this function can only be used to calculate
        // the objective function value but not the derivatives
        arma::vec Xglo, channelVec;
        //arma_mexprintf("channelNum",parent1.channelNum);
        //arma_mexprintf("channelLocNodes",parent1.channelLocNodes);
        for (std::size_t i = 0; i < parent1.channelNum.n_elem; i++)
        {
            if (fabs(channels.mcf(parent1.channelNum(i))) < MCFTOL)
                continue;
            channelVec = Xel.col(parent1.channelLocNodes(1,i))-Xel.col(parent1.channelLocNodes(0,i));
            for (std::size_t j = 0; j < gauss1.line.n_cols; j++) 
            { 
                Xglo = Xel.col(parent1.channelLocNodes(0,i)) + gauss1.line(0,j)*channelVec;
                locCoord = local_coord_2D(Xglo,Xel);
                shape_function_2D(N,DN,locCoord);
                uhloffset = arma::dot(N,UURel) - objOpt.offset;
                factor = 0.5*gauss1.line(1,j)*arma::norm(channelVec);
                objValEl += pow(uhloffset,objOpt.normp)*factor;
            } 
        }
        if (objOpt.calcOneNorm) 
        {
            for (std::size_t i = 0; i < gauss1.elem.n_cols; i++)
            {      
                shape_function_2D(N,DN,gauss1.elem(arma::span(0,gauss1.elem.n_rows-2),i));
                factor = detJ*gauss1.elem(gauss1.elem.n_rows-1,i);
                
                uh = arma::dot(N,UURel); 
                
                oneNormValEl += uh*factor;
                if (calcGrad) Fadj1normEl += N*factor;
            }
        }
    }
}   


}

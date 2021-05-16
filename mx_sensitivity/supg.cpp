/*
Created by Marcus Tan on 11/1/2014
Modified on 1/6/2015
Copyright 2014 University of Illinois 
Purpose: this function calculates the body source at a give position
INPUT:
    B: [dN1/dx,dN2/dx,...dNn/dx;
        dN1/dy,dN2/dy,...dNn/dy]';
        n is the number of shape functions
*/
#include "sensitivity.h"
#include "armadillo"
#include <math.h>
#include <mex.h>

namespace igfem
{
void streamwise_elem_length(double& he, 
                            arma::vec& Bsw,
                            const arma::vec& channelVec,
                            const arma::mat& B)
{
    Bsw = B*channelVec;
    
     double denom = arma::sum(arma::abs(Bsw));
     if (fabs(denom) < DENOMTOL)
     {
        //std::cerr << "streamwise_elem_length: denom close to zero" << std::endl;
        mexErrMsgIdAndTxt("streamwise_elem_length:denom_vanish", 
                           "streamwise_elem_length vanishing denominator");   
        he = 0.0;
        return;
     }
     he = 1.0/denom; // Note: this is half the streamwise length presented in papers
    
    return;
}

void supg_weighing_function(arma::vec& W,
                            const arma::mat& channelUnitVecs,
                            const arma::vec& N,
                            const arma::mat& B,
                            const arma::vec& mcf)
{
    arma::vec Bsw;
    W = N;
    double he;
    for (std::size_t i = 0; i < channelUnitVecs.n_cols;i++)
    {
        if (fabs(mcf(i)) < MCFTOL)
            continue;
        streamwise_elem_length(he,Bsw,sgn(mcf(i))*channelUnitVecs.col(i),B);
        W += he*Bsw;
    }
}

void D2_supg_weighing_function(arma::vec& D2W,
                               const arma::mat& channelUnitVecs,
                               const arma::mat& D2t,
                               const arma::uvec& chLocPaEnrichNodes,
                               const arma::vec& D2N,
                               const arma::mat& B,
                               const arma::mat& D2BchEn,
                               const arma::vec& mcf)
{
    D2W = D2N;
    arma::vec Bsw, unitVec, D2unitVec;
    double he,D2he;
    for (std::size_t i = 0; i < channelUnitVecs.n_cols;i++)
    {     
          if (fabs(mcf(i)) < MCFTOL)
            continue;
          unitVec = sgn(mcf(i))*channelUnitVecs.col(i);
          streamwise_elem_length(he,Bsw,unitVec,B);
          Bsw = arma::sign(Bsw);
          D2W(chLocPaEnrichNodes) += he*D2BchEn*unitVec;
          D2unitVec = sgn(mcf(i))*D2t.col(i);
          D2W += he*B*D2unitVec;
          D2he = -he*he*(arma::sum(Bsw(chLocPaEnrichNodes)
                           % (D2BchEn*unitVec)) 
                  +arma::sum(Bsw % (B*D2unitVec)));
          D2W += D2he*B*unitVec;
    }
}

}

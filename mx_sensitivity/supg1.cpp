/*
Created by Marcus Tan on 11/1/2014
Modified on 1/6/2015
Copyright 2014 University of Illinois 
Purpose: this function calculates the body source at a give position
INPUT:
    B: [dN1/dx,dN2/dx,...dNn/dx;
        dN1/dy,dN2/dy,...dNn/dy];
        n is the number of shape functions
*/
#include "sensitivity.h"
#include "armadillo"
#include <math.h>
#include <iostream>

namespace igfem
{
void streamwise_elem_length(double& he, 
                            arma::vec& Bsw,
                            const arma::vec& channelVec,
                            const arma::mat& B)
{
    Bsw = B*channelVec;
    //double denom = arma::sum(arma::abs(Bsw));
    //if (fabs(denom) < DENOMTOL)
    //    std::cerr << "streamwise_elem_length: denom close to zero" << std::endl;
    he = 1.0/arma::sum(arma::abs(Bsw));
    
    return;
}

void supg_weighing_function(arma::vec& W,
                            const arma::mat& channelUnitVecs,
                            const arma::vec& N,
                            const arma::mat& B)
{
    arma::vec Bsw;
    W = N;
    double he;
    for (std::size_t i = 0; i < channelUnitVecs.n_cols;i++)
    {
        streamwise_elem_length(he,Bsw,channelUnitVecs.col(i),B);
        W += he*Bsw;
    }
}

void D2_supg_weighing_function(arma::vec& D2W,
                               const arma::mat& channelUnitVecs,
                               const arma::mat& D2t,
                               const arma::uvec& chLocPaEnrichNodes,
                               const arma::vec& D2N,
                               const arma::mat& B,
                               const arma::mat& D2BchEn)
{
    D2W = D2N;
    arma::vec Bsw;
    double he,D2he;
    for (std::size_t i = 0; i < channelUnitVecs.n_cols;i++)
    {
          streamwise_elem_length(he,Bsw,channelUnitVecs.col(i),B);
          Bsw = arma::sign(Bsw);
          D2W(chLocPaEnrichNodes) += he*D2BchEn*channelUnitVecs.col(i);
          D2W += he*B*D2t.col(i);
          D2he = -he*he*(arma::sum(Bsw(chLocPaEnrichNodes)
                           % (D2BchEn*channelUnitVecs.col(i))) 
                  +arma::sum(Bsw % (B*D2t.col(i))));
          D2W += D2he*B*channelUnitVecs.col(i);
    }
}

}

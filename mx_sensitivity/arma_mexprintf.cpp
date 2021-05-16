/*
Created by Marcus Tan on 1/30/2016
Copyright 2016 University of Illinois 
Purpose: this function use mexPrintf to print Armadillo vectors and matrices since the Armadillo member function print
may not always work in Matlab
*/

#include "sensitivity.h"
#include <cstddef> // NULL, std::#include "mex.h"
#include "armadillo"
#include "mex.h"

namespace igfem
{
void arma_mexprintf(const char prefix[], const arma::vec& V)
{
    mexPrintf("%s\n",prefix);
    for (std::size_t i = 0; i < V.n_elem; i++)
        mexPrintf("%8.4g\n",V(i));
    mexPrintf("\n");
    return;
}

void arma_mexprintf(const char prefix[], const arma::uvec& V)
{
    mexPrintf("%s\n",prefix);
    for (std::size_t i = 0; i < V.n_elem; i++)
        mexPrintf("%i\n",V(i));
    mexPrintf("\n");
    return;
}
void arma_mexprintf(const char prefix[], const arma::mat& V)
{
    mexPrintf("%s",prefix);
    for (std::size_t i = 0; i < V.n_rows; i++)
    {
        mexPrintf("\n");
        for (std::size_t j = 0; j < V.n_cols; j++)
        {
            mexPrintf("%8.4g ",V(i,j));
        }
    }
    mexPrintf("\n");
    return;
}

void arma_mexprintf(const char prefix[], const arma::umat& V)
{
    mexPrintf("%s",prefix);
    for (std::size_t i = 0; i < V.n_rows; i++)
    {
        mexPrintf("\n");
        for (std::size_t j = 0; j < V.n_cols; j++)
        {
            mexPrintf("%i ",V(i,j));
        }
    }
    mexPrintf("\n");
    return;
}

}


/*
Created by Marcus Tan on 8/17/2014
Modified on 12/28/2014
Copyright 2014 University of Illinois 
Purpose: this function calculates shape functions, B matrix of the child 
        and parent elements as well as the Jacobian of the child element
        mapping to the global space
INPUT:
    paLocOrigNodes: parent local number of original nodes
    chLocPaEnrichNodes: child local number of enrichnment nodes wrt parent
    chLocEnrichNodes: child local number of enrichment nodes wrt itself
    shape: TRIANGULAR or QUADRILATERAL
    nNodes: total number of nodes of IGFEM element
    Xel: dim x nOriginalNodes matrix of original element coordinates
    Xch: dim x nChildNodes of child element coordinates
    locCoord: a vector of the local coordinates
OUTPUT: 
    N: shape functions
    B: derivative of shape functions wrt global coordinates 
        (note: the size of the matrix is nNodesPerElem x dim
                and it is the transpose of that in mx_FEM)
    invJ: inverse of Jacobian
    detJ: determinant of Jacobian
INPUT_OUTPUT: 
    DNch: nChildNodes x dim DN matrix of child element 
    Bch: nChildNodes of child element x dim B matrix of child element
    Bel: B matrix of original element
    calcBJch: flag indicating whether the B matrix of the child element should be calculated. 
              if true, it calculates the matrix and returns a false if child element is triangular
    calcBel: flag indicating whether the B matrix of the parent element should be calculated. 
             if true, it calculates the matrix and returns a false
*/
#include "sensitivity.h"
#include "armadillo"
#include <cstddef>

namespace igfem
{
void child_element_NDNBJ(arma::vec& N,
                       arma::mat& B,
                       arma::mat& DNch,
                       arma::mat& Bch,
                       arma::mat& DNel,
                       arma::mat& Bel,
                       arma::mat& invJ,
                       double& detJ,
                       bool& calcBJch,
                       bool& calcBel,
                       arma::vec& Xglo,
                       std::size_t nNodes,
                       const arma::uvec& paLocOrigNodes,
                       const arma::uvec& chLocPaEnrichNodes, 
                       const arma::uvec& chLocEnrichNodes, 
                       const arma::mat& Xel,
                       const arma::mat& Xch,
                       const arma::vec& locCoord,
                       shapes shape)
{
   
    N = arma::zeros<arma::vec>(nNodes);
    B = arma::zeros<arma::mat>(nNodes,Xel.n_rows);

    // calculate shape functions, B matrix and Jacobian of child element
    arma::vec Nch; // child shape function
    shape_function_2D(Nch,DNch,locCoord,shape);
    if (calcBJch)
    {
        arma::mat Jch = Xch*DNch;
        detJ = arma::det(Jch);
     
        if (detJ > JACTOL)
        {
            invJ = arma::inv(Jch);
            Bch = DNch*invJ;
        }
        else
        {
            invJ = arma::zeros<arma::mat>(Xel.n_rows,Xel.n_rows);
            Bch = arma::zeros<arma::mat>(DNch.n_rows,DNch.n_cols);
        }
        if (shape == TRIANGLE)     
            calcBJch = false;
    }
    N(chLocPaEnrichNodes) = Nch(chLocEnrichNodes);
    B.rows(chLocPaEnrichNodes) = Bch.rows(chLocEnrichNodes);
    
    // calculate shape functions, B matrix of nodes belonging to parent element
    Xglo = Xch*Nch;
    arma::vec locCoordPa = local_coord_2D(Xglo,Xel);
  
    arma::vec Nel;
    shape_function_2D(Nel,DNel,locCoordPa);
    
    if (calcBel)
    {
        arma::mat Jel = Xel*DNel;
        Bel = arma::trans(arma::solve(Jel.t(),DNel.t()));
        calcBel = false;
    }
   
    N(paLocOrigNodes) = Nel;
    B.rows(paLocOrigNodes) = Bel;
    
}
}

/*
Created by Marcus Tan on 8/7/2014
Update on 10/25/2014
Copyright 2014 University of Illinois 
Purpose: header file for the assembly of stiffness matrix
*/
#ifndef SENSITIVITY_H_
#define SENSITIVITY_H_
#include <cstddef> // NULL, std::size_t
#include "armadillo"

#define JACTOL 1e-13 // tolerance below which Jacobian is considered zero
#define N_ORIGINAL_NODES 3 // number of original nodes is 3, assuming triangular element is used 

//#define SHOW_DK_DP

namespace igfem
{

enum ObjTypes {PNORMTEMP,AVETEMP,COMPLIANCE,OBJLAST};
const int nObjTypes = PNORMTEMP - OBJLAST;
struct ObjOptions
{
    //ObjOptions(): normpNormalization(1.0){};
    ObjTypes type;
    double normp;
    //double normpNormalization;
    std::size_t nDesignParams;
};


// types of elements
enum shapes {TRIANGLE,QUADRILATERAL};
struct child
{
    shapes shape; // child element type
    arma::uvec locPaNodes; //local numbers of child element nodes relative to parent nodes
    arma::uvec locPaEnNodes; //local numbers of child element enrichment nodes rel to parent nodes
    arma::uvec locEnNodes; //local numbers of child element enrichment nodes rel to itself
    arma::uvec channelNum;
    arma::umat channelNodes;
    arma::umat channelLocNodes;
    // arma::mat channelNurbsParam;
};

enum parentType {REGULAR,IGFEM};
struct parent
{
    parent() : type(REGULAR){};
    parentType type;
    arma::uvec nodes; 
    arma::field<child> children;
    double conductivity;
    // designParamNum: a vector of design parameter numbers that have 
    // non-zero velocties in any child elem
    // corresponding to each slice of vel
    arma::uvec designParamNum;
    // vel: a ndim x nEnrichmentNodes x designParamNum.n_elem cube matrix
    arma::cube vel;    
};

enum modelType{MEAN_TEMP,CONST_HEAT};
struct chanNetwork
{
    modelType model;
    arma::vec mcf; // vector of mass flow rate x heat capacity of all channels
    // following variables are only used for the constant heat flux model
    //double Tin; // channel inlet temperature
    //arma::vec kapf; // thermal conductivity of fluid
    //arma::vec eigvalsq;
    //arma::vec CR1s;
    //arma::vec lengths; // lengths of all channels
    arma::mat D2mcf; // derivative of (mass flow rate x heat capacity) 
                        // wrt to each design parameter
};


struct gauss
{
    arma::mat elem; // a (number of dimensions+1)x number of gauss points matrix for integration over elements
    arma::mat quadElem;
    arma::mat line;  // a 2 x number of gauss points matrix for integration over lines. 
};

struct Neumann
{
    Neumann()
    {
        surf = -1;
    }
    int surf; // default not a Neumann element
    double val;
};

struct Convection
{
        double coef;
        double Tref;
};

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
                                    bool supg); //
                                    //const arma::vec& UP

void regular_element_pseudo_adjoint_forces(arma::mat& FpseudoEl,
                                           arma::vec& FadjEl,
                                           arma::vec& gradObjRemEl,
                                           double& objVal,
                                           const ObjOptions& objOpt,
                                           const arma::mat& Xe,
                                           double elemHeatSource,
                                           const gauss& gauss1,
                                           const arma::vec& UURel);

template <typename T> 
inline int sgn(T val) {
        return (T(0) < val) - (val < T(0));
}

#ifdef SHOW_DK_DP
void IGFEM_element_pseudo_adjoint_forces(arma::mat& FpseudoEl,
                                         arma::vec& FadjEl,
                                         arma::vec& gradObjRemEl,
                                         double& objVal,
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
                                         bool supg);
                                         //const arma::vec& UPel	
#else
void IGFEM_element_pseudo_adjoint_forces(arma::mat& FpseudoEl,
                                         arma::vec& FadjEl,
                                         arma::vec& gradObjRemEl,
                                         double& objVal,
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
                                         bool supg);
  
#endif
void shape_function_2D(arma::vec& N, arma::mat& DN, const arma::vec& locCoord); // triangular shape function

void shape_function_2D(arma::vec& N, arma::mat& DN, const arma::vec& locCoord, 
                       shapes shape); // triangular or quadrilateral shape functions


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
                       shapes shape);



arma::vec local_coord_2D(const arma::vec& X, const arma::mat& Xel);

arma::vec local_coord_2D_along_edge(const arma::vec& X, const arma::mat& Xel, shapes shape, int edge);

int edge_number(const arma::uvec& twoNodes, shapes shape);

// for supg calculation
void streamwise_elem_length(double& he, 
                            arma::vec&Bsw1, 
                            const arma::vec& channelVec,
                            const arma::mat& B);
							
void supg_weighing_function(arma::vec& W,
                            const arma::mat& channelUnitVecs, 
                            const arma::vec& N,
                            const arma::mat& B,
                            const arma::vec& mcf);

void D2_supg_weighing_function(arma::vec& D2W,
                               const arma::mat& channelUnitVecs,
                               const arma::mat& D2t,
                               const arma::uvec& chLocPaEnrichNodes,
                               const arma::vec& D2N,
                               const arma::mat& B,
                               const arma::mat& D2BchEn,
                               const arma::vec& mcf);

void body_source(double& fb,arma::vec& Dfb,const arma::vec& X);

}
#endif

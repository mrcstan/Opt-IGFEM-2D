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
#define DENOMTOL 1e-13 // tolerance below which a denominator is considered zero
#define N_ORIGINAL_NODES 3 // number of original nodes is 3, assuming triangular element is used 
#define MCFTOL 1e-13 // tolerance below which the mass flow rate x heat capacity is considered zero 
                    // (should be consistent with that in mx_FEM/assemble.h

//#define SHOW_DK_DP // directive to show derivative of stiffness matrix and load vector wrt design parameters
// NOTE: SHOW_DK_DP directive must enabled simultaneously in the following source codes:
//       assemble_pseudo_adjoint_forces.cpp
//       IGFEM_element_pseudo_adjoint_forces.cpp

namespace igfem
{

//enum ObjTypes {PNORMTEMP,AVETEMP,COMPLIANCE,OBJLAST};
//const int nObjTypes = PNORMTEMP - OBJLAST;

enum IntDomainTypes {WHOLE,CHANNEL};
struct ObjOptions
{
    //ObjOptions(): normpNormalization(1.0){};
    //ObjTypes type;
    //double normpNormalization;
    double normp; // p value in the p-norm
    double offset; // subtract offset from temperature before calculating p-norm
    std::size_t nDesignParams; // number of design parameters
    bool calcOneNorm; // a flag that indicates whether besides the objective function, 
                     // the 1-norm should also be calculated
    IntDomainTypes intDomain; // integration domain (can either be the whole 2D domain or along the channels)
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
    arma::umat channelLocNodes;
    arma::uvec channelNum;
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

enum tempUnit {KELVIN,CELSIUS};
struct Convection
{
        double coef;
        double Tamb;
        double epsSB;
        tempUnit Tunit;
        bool linearRad;
};

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
                                    bool calcGrad); 

void regular_element_adjoint_forces(arma::vec& FadjEl,
                                    double& objVal,
                                    arma::vec& Fadj1normEl,
                                    double& oneNormVal,
                                    const ObjOptions& objOpt,
                                    const arma::mat& Xe,
                                    double elemHeatSource,
                                    const gauss& gauss1,
                                    const parent& parent1,
                                    const chanNetwork& channels,
                                    const arma::vec& UURel,
                                    bool calcGrad);

template <typename T> 
inline int sgn(T val) {
        return (T(0) < val) - (val < T(0));
}

#ifdef SHOW_DK_DP
void IGFEM_element_pseudo_adjoint_forces(arma::mat& FpseudoEl,
                                         arma::vec& FadjEl,
                                         arma::vec& gradObjRemEl,
                                         double& objVal,   
                                         arma::vec& Fadj1normEl,
                                         arma::vec& grad1normRemEl,
                                         double& oneNormVal,
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
                                         bool calcGrad);
#else
void IGFEM_element_pseudo_adjoint_forces(arma::mat& FpseudoEl,
                                         arma::vec& FadjEl,
                                         arma::vec& gradObjRemEl,
                                         double& objVal,
                                         arma::vec& Fadj1normEl,
                                         arma::vec& grad1normRemEl,
                                         double& oneNormVal,
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
                                         bool calcGrad);
  
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

void arma_mexprintf(const char prefix[], const arma::vec& V);
void arma_mexprintf(const char prefix[], const arma::uvec& V);
void arma_mexprintf(const char prefix[], const arma::mat& V);
void arma_mexprintf(const char prefix[], const arma::umat& V);

}
#endif

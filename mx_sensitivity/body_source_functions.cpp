/*
Created by Marcus Tan on 2/12/2015
Copyright 2014 University of Illinois 
Purpose: this function calculates the body source at a give position
*/
#include "sensitivity.h"
#include "armadillo"
#include <iostream>
#include <math.h>

#define X0Y0 0.2489
#define X1Y0 23.42
#define X0Y1 21.87
#define X2Y0 -485.7
#define X1Y1 -16.96
#define X0Y2 -356.8
#define X3Y0 4317
#define X2Y1 170.4
#define X1Y2 276.7
#define X0Y3 2422
#define X4Y0 -14060
#define X3Y1 -791
#define X2Y2 -490
#define X1Y3 -972.2
#define X0Y4 -5849

// #define FOR_VALIDATION
// #define DOUBLE_LOCALIZED_HEAT
// #define GRIDLIKE_CHANNELS
#define GRIDLIKE_NETWORK_LOCALIZED_HEAT

// To use any of the distributed heat sources in this source code,
// the directive SOURCE_FUNCTION must be uncommented
// in the following source codes:
// IGFEM_element_pseudo_adjoint_forces.cpp
namespace igfem
{
#ifdef FOR_VALIDATION
void body_source(double& fb,arma::vec& Dfb,const arma::vec& X)
{
    Dfb.set_size(2);       
    // Stephen's experiment, non-uniform heat source
    fb = 500*(X0Y0 + X1Y0*X(0) + X0Y1*X(1) + X2Y0*X(0)*X(0) + X1Y1*X(0)*X(1) + X0Y2*X(1)*X(1)
             + X3Y0*pow(X(0),3.0) + X2Y1*X(0)*X(0)*X(1) + X1Y2*X(0)*X(1)*X(1) + X0Y3*pow(X(1),3.0)
             + X4Y0*pow(X(0),4.0) + X3Y1*pow(X(0),3.0)*X(1) + X2Y2*X(0)*X(0)*X(1)*X(1)
             + X1Y3*X(0)*pow(X(1),3.0) + X0Y4*pow(X(1),4.0));
    Dfb(0) = 500*(X1Y0 + 2*X2Y0*X(0) + X1Y1*X(1) 
                + 3*X3Y0*X(0)*X(0) + 2*X2Y1*X(0)*X(1) + X1Y2*X(1)*X(1)
                + 4*X4Y0*pow(X(0),3.0) + 3*X3Y1*X(0)*X(0)*X(1) + 2*X2Y2*X(0)*X(1)*X(1)
                + X1Y3*pow(X(1),3.0));
    Dfb(1) = 500*(X0Y1 + X1Y1*X(0) + 2*X0Y2*X(1)
                 + X2Y1*X(0)*X(0) + 2*X1Y2*X(0)*X(1) + 3*X0Y3*X(1)*X(1)
                 + X3Y1*pow(X(0),3.0) + 2*X2Y2*X(0)*X(0)*X(1)
                 + 3*X1Y3*X(0)*X(1)*X(1) + 4*X0Y4*pow(X(1),3.0));
}
#endif

#ifdef DOUBLE_LOCALIZED_HEAT
void body_source(double& fb,arma::vec& Dfb,const arma::vec& X)
{ 
    const double ro = 0.015;
    const double xo = 0.04;
    const double yo = 0.04;
    const double Qo = 250*0.15*0.2*pow(15/(16*ro),2);
    const double xmin = xo - ro;
    const double xmax = xo + ro;
    const double ymin = yo - ro;
    const double ymax = yo + ro;
    
    fb = 0.0;
    Dfb.set_size(2);
    Dfb(0) = 0.0;
    Dfb(1) = 0.0;
    
    if (X(0) >= xmin && X(0) <= xmax && X(1) >= ymin && X(1) <= ymax)
    {
        fb = Qo*pow(1-pow((X(0)-xo)/ro,2),2)*pow(1-pow((X(1)-yo)/ro,2),2);
        Dfb(0) = 4*Qo/(ro*ro)*(xo-X(0))*(1-pow((X(0)-xo)/ro,2))*pow(1-pow((X(1)-yo)/ro,2),2);
        Dfb(1) = 4*Qo/(ro*ro)*(yo-X(1))*(1-pow((X(1)-yo)/ro,2))*pow(1-pow((X(0)-xo)/ro,2),2);
    }
    
    const double r1 = 0.015;
    const double x1 = 0.11;
    const double y1 = 0.16;
    const double Q1 = 250*0.15*0.2*pow(15/(16*r1),2);
    const double xmin1 = x1 - r1;
    const double xmax1 = x1 + r1;
    const double ymin1 = y1 - r1;
    const double ymax1 = y1 + r1;
    
    if (X(0) >= xmin1 && X(0) <= xmax1 && X(1) >= ymin1 && X(1) <= ymax1)
    {
        fb += Q1*pow(1-pow((X(0)-x1)/r1,2),2)*pow(1-pow((X(1)-y1)/r1,2),2);
		Dfb(0) += 4*Q1/(r1*r1)*(x1-X(0))*(1-pow((X(0)-x1)/r1,2))*pow(1-pow((X(1)-y1)/r1,2),2);
        Dfb(1) += 4*Q1/(r1*r1)*(y1-X(1))*(1-pow((X(1)-y1)/r1,2))*pow(1-pow((X(0)-x1)/r1,2),2);
    }
}
#endif


#ifdef GRIDLIKE_CHANNELS
void body_source(double& fb,arma::vec& Dfb,const arma::vec& X)
{    
    // For blockage tolerant/redundancy project
    fb = 0.0;
    Dfb.set_size(2);
    Dfb(0) = 0.0;
    Dfb(1) = 0.0;
    
    const double qave = 2000;
    fb = qave*(0.188625 + 127*X(0) + 99.54*X(1) 
               - 7583*X(0)*X(0) - 4038*X(0)*X(1) - 5169*X(1)*X(1)
              + 1.919E5*pow(X(0),3.0) + 2.205E5*X(0)*X(0)*X(1) + 6.490E4*X(0)*X(1)*X(1) 
              + 1.228E5*pow(X(1),3.0) - 1.700E6*pow(X(0),4.0) - 5.539E6*pow(X(0),3.0)*X(1) 
              - 1.921E6*X(0)*X(0)*X(1)*X(1) + 1.263E6*X(0)*pow(X(1),3.0) - 2.156E6*pow(X(1),4.0)
              - 2.348E6*pow(X(0),5.0) + 5.443E7*pow(X(0),4.0)*X(1) - 8.962E6*pow(X(0),3.0)*pow(X(1),2.0)
              + 4.017E7*pow(X(0),2.0)*pow(X(1),3.0) - 4.462E7*X(0)*pow(X(1),4.0) + 1.982E7*pow(X(1),5.0));
    Dfb(0) = qave*(         127
               - 2*7583*X(0) - 4038*X(1) 
              + 3*1.919E5*pow(X(0),2.0) + 2*2.205E5*X(0)*X(1) + 6.490E4*X(1)*X(1) 
                                        - 4*1.700E6*pow(X(0),3.0) - 3*5.539E6*pow(X(0),2.0)*X(1) 
              - 2*1.921E6*X(0)*X(1)*X(1) + 1.263E6*pow(X(1),3.0)
              - 5*2.348E6*pow(X(0),4.0) + 4*5.443E7*pow(X(0),3.0)*X(1) - 3*8.962E6*pow(X(0),2.0)*pow(X(1),2.0)
              + 2*4.017E7*X(0)*pow(X(1),3.0) - 4.462E7*pow(X(1),4.0));
    Dfb(1) = qave*(              99.54 
                               - 4038*X(0) - 2*5169*X(1)
                               + 2.205E5*X(0)*X(0) + 2*6.490E4*X(0)*X(1)
              + 3*1.228E5*pow(X(1),2.0)              - 5.539E6*pow(X(0),3.0) 
              - 2*1.921E6*X(0)*X(0)*X(1) + 3*1.263E6*X(0)*pow(X(1),2.0) - 4*2.156E6*pow(X(1),3.0)
                                         + 5.443E7*pow(X(0),4.0) - 2*8.962E6*pow(X(0),3.0)*X(1)
              + 3*4.017E7*pow(X(0),2.0)*pow(X(1),2.0) - 4*4.462E7*X(0)*pow(X(1),3.0) + 5*1.982E7*pow(X(1),4.0)); 
    
}
#endif

#ifdef GRIDLIKE_NETWORK_LOCALIZED_HEAT
void body_source(double& fb,arma::vec& Dfb,const arma::vec& X)
{
    // For blockage tolerant/redundancy project
    const double ro = 0.0125;
    const double xo = 0.05;
    const double yo = 0.05;
    const double Qo = 2000*0.05*0.05*pow(15/(16*ro),2);
    const double xmin = xo - ro;
    const double xmax = xo + ro;
    const double ymin = yo - ro;
    const double ymax = yo + ro;
    
    fb = 0.0;
    Dfb.set_size(2);
    Dfb(0) = 0.0;
    Dfb(1) = 0.0;
    
    if (X(0) >= xmin && X(0) <= xmax && X(1) >= ymin && X(1) <= ymax)
    {
        fb = Qo*pow(1-pow((X(0)-xo)/ro,2),2)*pow(1-pow((X(1)-yo)/ro,2),2);
        Dfb(0) = 4*Qo/(ro*ro)*(xo-X(0))*(1-pow((X(0)-xo)/ro,2))*pow(1-pow((X(1)-yo)/ro,2),2);
        Dfb(1) = 4*Qo/(ro*ro)*(yo-X(1))*(1-pow((X(1)-yo)/ro,2))*pow(1-pow((X(0)-xo)/ro,2),2);
    }
 
}
#endif
}

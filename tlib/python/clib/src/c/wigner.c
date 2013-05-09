//==================================================================================================
//
// Purpose: compute Wigner 3J symbols for special cases that are needed for aberration problem.
//
//==================================================================================================
//
// Author: Jens Chluba (Johns Hopkins University)
//
// first implementation: Nov  2012
// last modification   : Nov  2012
//
//==================================================================================================

//#include <iostream>
//#include <string.h>
//#include <fstream.h>
#include <math.h>
//#include <vector>

#include "utilis.h"
#include "wigner.h"

//using namespace std;

//==================================================================================================
//
// Wigner 3J symbol for m1=m2=m3=0. This symbol is needed for fixing the normalization of the
// recursion relations
//
//==================================================================================================
double ThreeJ_l1l2l3(int l1, int l2, int l3)
{
    if(fabs(l1-l2)>l3 || l1+l2<l3) return 0;
    
    int J=l1+l2+l3, g=J/2;
    double phase= g%2 ? -1 : 1;
    
    if(J%2) return 0;
    
    double a=log10factorial(g)-log10factorial(g-l1)
            -log10factorial(g-l2)-log10factorial(g-l3);
    
    double b=log10factorial(J-2*l1)+log10factorial(J-2*l2)
            +log10factorial(J-2*l3)-log10factorial(J+1);
    
    return phase*pow(10.0, 0.5*b+a);
}

/*
//==================================================================================================
// C and D functions of Luscombe & Luban,  1998
//==================================================================================================
double C_func(int m, int j1, int j2, int j3, int m1)
{
    return sqrt( 1.0*(j2-m+1)*(j2+m)*(j3-m1-m+1)*(j3+m1+m) );
}

double D_func(int m, int j1, int j2, int j3, int m1)
{
    return (j2+1)*j2+(j3+1)*j3-(j1+1)*j1-2*m*(m+m1);
}

//==================================================================================================
//
// Wigner 3J symbol for m1=-m2 and m3=0. 
//
//==================================================================================================
void Wigner_3J_l1l2l3_m_mm(int l1, int l2, int l3, 
                           vector<double> &results, 
                           int &mmin, int &mmax)
{
    results.clear();
    
    if(l1<0 || l2<0 || l3<0)
    { 
        cerr << " Something is wrong with l1, l2, or l3 " << endl; 
        exit(0); 
    }
    
    // check if l1, l2, l3 follow triangular relation
    if(fabs(l1-l2)>l3 || l1+l2<l3)
    {
        mmin=mmax=0;
        results.push_back(0.0);
        
        return;
    }
        
    // for the recursions m1=0, l1==j2, l2==j3 and l3==j1
    int j1=l3;
    int j2=l1;
    int j3=l2;
    
    //===============================================================================
    // memory and m-range
    //===============================================================================
    mmin=max(-j2,-j3);
    mmax=min(j2,j3);
    int mtot=mmax-mmin+1;
    
    double W0=ThreeJ_l1l2l3(l1, l2, l3);
    
    results.resize(mmax-mmin+1);
 
    //===============================================================================
    // simplest case
    //===============================================================================
    if(mmax==0)
    {
        results[0]=W0;
        return;
    }
    
    //===============================================================================
    // start recursions from mmin for m=mmin --> 0
    //===============================================================================
    results[0]=1.0;

    double X=C_func(mmin+1, j1, j2, j3, 0);
    double Y=D_func(mmin  , j1, j2, j3, 0);

    results[1]=-results[0]*Y/X; 

    double Z=X;
    
    for(int m=mmin+1, ind=2; m<1; m++, ind++)
    {
        X=C_func(m+1, j1, j2, j3, 0);
        Y=D_func(m  , j1, j2, j3, 0);
        
        results[ind]=-(results[ind-1]*Y+results[ind-2]*Z)/X;
        
        Z=X;
    }
     
    //===============================================================================
    // start recursions from mmax for m=mmax --> 0
    //===============================================================================
    results[mtot-1]=1.0;
    
    Z=C_func(mmax, j1, j2, j3, 0);
    Y=D_func(mmax, j1, j2, j3, 0);
    
    results[mtot-2]=-results[mtot-1]*Y/Z;
       
    X=Z;
    
    for(int m=mmax-1, ind=mtot-3; m>1; m--, ind--)
    {
        Z=C_func(m, j1, j2, j3, 0);
        Y=D_func(m, j1, j2, j3, 0);
        
        results[ind]=-(results[ind+2]*X+results[ind+1]*Y)/Z;
        
        X=Z;
    }

    double fac=W0/results[-mmin], phase=1.0;

    if(W0==0)
    {
        fac=0.0;
        
        for(int k=0; k<-mmin; k++) 
            fac+=(results[k]*results[k]+results[mtot-1-k]*results[mtot-1-k]);
            
        fac*=2.0*j1+1.0;
        fac=sqrt(1.0/fac)*onetothek(j2-j3);
        
        X=C_func(1, j1, j2, j3, 0);
        Z=C_func(0, j1, j2, j3, 0);
        
        phase= X*results[-mmin+1]*Z*results[-mmin-1] > 0.0 ? -1 : 1;
    }

    for(int k=0; k<-mmin; k++)
    {
        results[k]       *=fac;
        results[mtot-1-k]*=fac*phase;
    }
    
    results[-mmin]=W0;
        
    
    return;
}

//==================================================================================================
//==================================================================================================

*/

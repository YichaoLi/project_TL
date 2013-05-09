""">>> lg_angle_2pt.py <<<

This module implements the wide-angle redshift formulae
from Papai & Szapudi, arXiv:0802.2940.  Written by Peter Papai.

Current revision:
    ID:         $ID: $
    Date:       $Date: 2008-01-07 14:53:36 -1000 (Mon, 07 Jan 2008) $
    Revision:   $Revision: 56 $

(C) 2008 The CosmoPy Team (see Copyright for details)
"""

import hankel
import utils
import pylab as M
import sys


def gg1(beta,alpha):
    th2 = (alpha+2.*beta)/2.
    return M.sin(th2)/M.sin(alpha)


def gg2(beta,alpha):
    th1 = (2.*beta-alpha)/2.
    return M.sin(th1)/M.sin(alpha)


# fpower is a two column array with the 'k' the first and 'P(k)' in the second column;
# alpha is the array of opening angles; 
# beta is the array of angles between the line of sight and the orientation of the pair;
# r is the array of separations;
# fout is the output file;
# f is f;
# the rest is the parameters for the Hankel transformation;

def lg_angle_corr (fpower, alpha, beta, r, fout, f=0.485,n=1000,h=0.001):
    k = fpower[:,0]
    pk = fpower[:,1]
    # interpolate/extrapolate P(k) (power law extrapolation)
    pki = utils.LogInterp(k,pk)
    pkii = utils.LogInterp(k,pk/k)
    pkiii = utils.LogInterp(k,pk/(k**2))
    # init hankel transform
    h3d0 = hankel.Hankel3D(0.5,n)
    h3d2 = hankel.Hankel3D(2.5,n)
    h3d4 = hankel.Hankel3D(4.5,n)
    h3d1 = hankel.Hankel3D(1.5,n)
    h3d3 = hankel.Hankel3D(3.5,n)
    # hankel transform 
    xi0 = h3d0.transform(pki.f,r,n,h)
    xi2 = h3d2.transform(pki.f,r,n,h)   
    xi4 = h3d4.transform(pki.f,r,n,h)
    xii1 = h3d1.transform(pkii.f,r,n,h)
    xii3 = h3d3.transform(pkii.f,r,n,h)
    xiii0 = h3d0.transform(pkiii.f,r,n,h)
    xiii2 = h3d2.transform(pkiii.f,r,n,h)
    # calculation of the 'a'-s and 'b'-s
    a0 = (1+2*f/3+2*f*f/15)*xi0-(f/3+2*f*f/21)*xi2+(3*f*f/140)*xi4
    a02 = -(f/2+3*f*f/14)*xi2+(f*f/28)*xi4
    a22 = (f*f/15)*xi0-(f*f/21)*xi2+(19*f*f/140)*xi4
    b22 = (f*f/15)*xi0-(f*f/21)*xi2-(4*f*f/35)*xi4
    a10 = ((2*f+4.*f*f/5.)*xii1/r-f*f*xii3/(5.*r))  # this is a10 tilde in the paper
    a11 = (4*f*f*xiii0/(3.*r*r)-8*f*f*xiii2/(3.*r*r))
    a21 = (3*f*f*xii3/(5.*r)-2*f*f*xii1/(5.*r))
    b11 = (4*f*f*xiii0/(3.*r*r)+4*f*f*xiii2/(3.*r*r))
    b21 = -(2*f*f*xii1/(5.*r)+2*f*f*xii3/(5.*r))
    # this is the output file	
    fi = open(fout,'w')
    for i in alpha:
         for j in beta:
	     g1 = gg1(j,i)
	     g2 = gg2(j,i)
             for k in range(0,M.size(r)):
                th1 = (2.*j-i)/2.
                th2 = (i+2.*j)/2.
                original = a0[k]+a02[k]*M.cos(2.*th1)+a02[k]*M.cos(2.*th2)+a22[k]*M.cos(2.*th1)*M.cos(2.*th2)+b22[k]*M.sin(2.*th1)*M.sin(2.*th2)
                new = a10[k]*(M.cos(th1)/g1-M.cos(th2)/g2)+a11[k]*M.cos(th1)*M.cos(th2)/(g1*g2)+a21[k]*(M.cos(2.*th1)*M.cos(th2)/g2-M.cos(th1)*M.cos(2.*th2)/g1)+b11[k]*M.sin(th1)*M.sin(th2)/(g1*g2)+b21[k]*(M.sin(2.*th1)*M.sin(th2)/g2-M.sin(th1)*M.sin(2.*th2)/g1)		
                if (2.*j>i): 
                        fi.write(str(i)+" "+str(j)+" "+str(r[k])+" "+str(new + original)+" "+"\n")
    fi.close()










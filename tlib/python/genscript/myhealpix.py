import numpy as np
import pylab as pl
import healpy as hp




def map_rotation(m, ang=None, coord=None, inv=None, deg=True, eulertype='X'):
    """ Return a rotated map with Euler angle `ang`
        `m`:   origin map
        `ang`: 3-list of Euler angle
	`eulertype`: 
    """

    nside=hp.get_nside(m)

    m_new=np.zeros(len(m) )
    pixcoord=np.arange(len(m))
    theta_after, phi_after = hp.pix2ang(nside, pixcoord)

    if coord!=None:
        rot=hp.Rotator(coord=coord)
    else:
        rot=hp.Rotator(rot=ang, inv=inv, deg=deg, eulertype=eulertype) 

    theta, phi=rot.I(theta_after, phi_after)
    m_new=hp.get_interp_val(m, theta, phi)


    return m_new





def map_mark(m, ang, deg=True, value=999, dnslice=8, nest=False):
    """ `m`:       original map
        `ang`:     position [theta, phi] of marked region
	`value`:   value of marked region
    """

    nslice=dnslice
    emp=np.zeros(hp.nside2npix(nslice))

    if deg==True:
        theta, phi= ang[0]*np.pi/180., ang[1]*np.pi/180.
    else:
        theta, phi=ang[0], ang[1]

    ind=hp.ang2pix(nslice, theta, phi)
    emp[ind]=value

    map_mark=hp.ud_grade(map_in=emp, nside_out=hp.get_nside(m))
    m_new=m+map_mark

    return m_new



def thetaphi2bl(theta, phi):
    ''' Return Galactic coordinate (l, b) '''

    return  


def bl2thetaphi(b, l, deg=True):

    if deg==True:
        return 90.-b, l
    else:
        raise Exception



def euler_angle_toZ(b, l):
    return l-90., 90.-b, 0.

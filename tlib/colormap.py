import numpy as np
import pylab as pl
from matplotlib import colors




def Colormap_Base_3color(color1, color2, scale1, scale2, cname='mycmap'):
    cmap=colors.LinearSegmentedColormap.from_list(cname,[color1, color2], 16)
    cmap._init()
    alphas=np.linspace(scale1, scale2, cmap.N+3)
    cmap._lut[:,-1] = alphas 

    return  cmap


def colormap_WGB(scale_1=0.1, scale_2=0.7):
    c1, c2='Green', 'Blue'

    cm = Colormap_Base_3color(c1, c2, scale_1, scale_2, cname='my_WGB')

    return cm




def colormap(color_scheme):
    sname_wgb=['White_Green_Blue', 'WGB']

    if color_scheme in sname_wgb:
        return colormap_WGB()



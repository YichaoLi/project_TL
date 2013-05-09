""" >>> plot.py <<<
    General Ploting Module
    Last modification: 2012-07-29 (01:03:54)
"""
import numpy as np
import pylab as pl
from matplotlib import colors, ticker

import myutlis as mut


""" Color map generator
"""
def gen_colormap(c=None):
    """ get the color map of mask """

    if c==None:
        color1=colors.colorConverter.to_rgba('grey')
        color2=colors.colorConverter.to_rgba('black')
    else:
        pass

    #print color1, color2

    my_cmap=colors.LinearSegmentedColormap.from_list('mycmap',[color1, color2], 16)
    my_cmap._init()
    alphas=np.linspace(0.3, 0.7, my_cmap.N+3)
    my_cmap._lut[:,-1] = alphas 

    return  my_cmap





''' >>>class 'myfig' contains all methods to make oridinary 1-D plots<<< '''
class fig(object):
    def __init__(self):
        pass
    
    def plot1d(self, data_x, data_y, output, xlim=0, ylim=0, scale=['log', 'log']):
	pl.clf()
        if(len(data_x.shape)==1):
	    pl.plot(data_x[:], data_y[:], 'k-', linewidth=1.3)
        else:
	    for i in range(data_x.shape[0]):
	        pl.plot(data_x[i,:], data_y[i,:], 'k-', linewidth=1.3)

        if((xlim==0) | (ylim==0)):
            xlim=[data_x.min(), data_x.max() ]
            ylim=[data_y.min(), data_y.max() ]

	pl.xscale(scale[0])
	pl.yscale(scale[1])
	pl.xlim(xlim)
	pl.ylim(ylim)

	if(output!=None):
	    pl.savefig(output)


    def cont2d(self, data, axis, output, scale=['log', 'log'], save=True):
        pl.clf()
	nl_f, nl_l = 5., 20

	dat = np.ma.masked_where(data<=0, data)

        if((scale[0]=='linear') & (scale[1]=='linear')):
            maxx, minn  = np.max(data), np.min(data)
            level = np.arange(minn, maxx, (maxx-minn)/nl_f )
            contf =pl.contourf(axis[0,:], axis[1,:], dat, level, \
	                     cmap=pl.cm.get_cmap('Blues'), shading='gouraud')
        else:
            contf =pl.contourf(axis[0,:], axis[1,:], dat, locator=ticker.LogLocator(), \
	                     cmap=pl.cm.get_cmap('Blues'), shading='gouraud')


        cont1=pl.contour(axis[0,:], axis[1,:], dat, 3, colors='k', linewidth=3, 
	                  alpha=0.5, hold = 'on')

	pl.colorbar(contf)

	""" different levels
        maxx, minn  = np.log10(np.max(dat) ), np.log10(1e-50)
        level = np.power(10,  np.arange(minn, maxx, (maxx-minn)/nl_f ) )
        contf =pl.contourf(self.xaxis[findex], self.yaxis[findex], dat, level, \
			 norm=colors.LogNorm() ,cmap=pl.cm.get_cmap('Blues'))
	# contour lines
        contl =pl.contour(self.xaxis[findex], self.yaxis[findex], dat, nl_l, \
		    	    colors='k', linewidth=3, alpha=0.5, hold = 'on')
	"""

        pl.xlim([axis[0,:].min(), axis[0,:].max()])
        pl.ylim([axis[1,:].min(), axis[1,:].max()])

	pl.xscale(scale[0])
	pl.yscale(scale[1])
	if(save==True):
	    pl.savefig(output)



    def cont2d_uneq(self, data, axis, output, scale=['log', 'log'], save=True, pcolors='Blues'):
        pl.clf()
	nl_f, nl_l = 5., 20

	dat = np.ma.masked_where(data<=0, data)

        if((scale[0]=='linear') & (scale[1]=='linear')):
            maxx, minn  = np.max(data), np.min(data)
            level = np.arange(minn, maxx, (maxx-minn)/nl_f )
            contf =pl.contourf(axis[0], axis[1], dat, level, \
	                     cmap=pl.cm.get_cmap(pcolors), shading='gouraud')
        else:
            contf =pl.contourf(axis[0], axis[1], dat, locator=ticker.LogLocator(), \
	                     cmap=pl.cm.get_cmap(pcolors), shading='gouraud')


        cont1=pl.contour(axis[0], axis[1], dat, 3, colors='k', linewidth=3, 
	                  alpha=0.5, hold = 'on')

	pl.colorbar(contf)

	""" different levels
        maxx, minn  = np.log10(np.max(dat) ), np.log10(1e-50)
        level = np.power(10,  np.arange(minn, maxx, (maxx-minn)/nl_f ) )
        contf =pl.contourf(self.xaxis[findex], self.yaxis[findex], dat, level, \
			 norm=colors.LogNorm() ,cmap=pl.cm.get_cmap('Blues'))
	# contour lines
        contl =pl.contour(self.xaxis[findex], self.yaxis[findex], dat, nl_l, \
		    	    colors='k', linewidth=3, alpha=0.5, hold = 'on')
	"""

        pl.xlim([axis[0].min(), axis[0].max()])
        pl.ylim([axis[1].min(), axis[1].max()])

	pl.xscale(scale[0])
	pl.yscale(scale[1])
	if(save==True):
	    pl.savefig(output)





    def mesh2d(self, data, axis, output=None, mask=None, datscale='log', axiscale=['log', 'log'], pcolors='Greys', maskcolors='blue', malpha=0.2):
        pl.clf()

	#ppldat = np.ma.masked_where(data<=0, data)
	pldat= data #ppldat

        """
	if mask!=None:
	    #pldat=np.ma.masked_where((mask<=1e-2)&(mask>=-1e-2), ppldat)
	    pldat=np.ma.masked_where((mask>1e-2)|(mask<-1e-2), ppldat)
	"""


	if(datscale=='log'):
	    cnorm=colors.LogNorm()
	elif(datscale=='linear'):
	    cnorm=colors.NoNorm()
	else:
	    raise Exception

        if pcolors!=None:
            # get the underlying map first
            cm=pl.pcolormesh(axis[0,:], axis[1,:], pldat, cmap=pl.cm.get_cmap(pcolors),
	                      norm=cnorm) 
	else:
            cm=pl.pcolormesh(axis[0,:], axis[1,:], pldat, norm=cnorm) 
	#pl.colorbar(cm)

        if mask!=None:
            # if do the mask

            #pl.pcolormesh(axis[0,:], axis[1,:], maskdata, alpha=0.3,  
	    #          cmap=pl.cm.get_cmap(maskcolors)) #, shading='gouraud')
            ##pl.pcolormesh(axis[0,:], axis[1,:], maskdata, alpha=0.4 ) 

	    if isinstance(maskcolors, list):
	        # multi-components of mask
	        for i in range(len(maskcolors)):
	            maskdata=np.ma.masked_where((mask-i-1>1e-2)|(mask-i-1<-1e-2),mask)
		    if isinstance(malpha, list):
		        mask_alpha=malpha[i]
		    else:
		        mask_alpha=malpha
	            if isinstance(maskcolors[i], str):
                        pl.contourf(axis[0,:], axis[1,:], maskdata, alpha=mask_alpha, 
	                         cmap=pl.cm.get_cmap(maskcolors[i]))
		    else:
                        pl.contourf(axis[0,:], axis[1,:], maskdata, alpha=mask_alpha, 
	                          cmap=maskcolors[i])

            else:
	        maskdata=np.ma.masked_where((mask<=1e-2)&(mask>=-1e-2), mask)

	        if isinstance(maskcolors, str):
                    pl.contourf(axis[0,:], axis[1,:], maskdata, alpha=malpha, 
	                    cmap=pl.cm.get_cmap(maskcolors))
	        else:
                    pl.contourf(axis[0,:], axis[1,:], maskdata, alpha=malpha, 
	                    cmap=maskcolors)

	pl.xscale(axiscale[0])
	pl.yscale(axiscale[1])

	if output!=None:
	    pl.savefig(output)

	return



    def mesh2d_mcolor_mask(self, data, axis, output=None, mask=None, datscale='log', 
           axiscale=['log', 'log'], pcolors='Greys', maskcolors=None):

        """       >>> generate 2D mesh plot <<<
	"""

        pl.clf()
	fig=pl.figure()
	ax=fig.add_subplot(111)

	pldat=data  
	
        # get the color norm
	if(datscale=='log'):
	    cnorm=colors.LogNorm()
	elif(datscale=='linear'):
	    cnorm=colors.NoNorm()
	else:
	    raise Exception


        color1=colors.colorConverter.to_rgba('white')
        color2=colors.colorConverter.to_rgba('blue')
        color3=colors.colorConverter.to_rgba('yellow')
        my_cmap0=colors.LinearSegmentedColormap.from_list('mycmap0',[color1, color1, color2, color2, color2, color3, color3], 512) 
        my_cmap0._init()


        if pcolors!=None:
            cm=ax.pcolormesh(axis[0,:], axis[1,:], pldat, cmap=pl.cm.get_cmap(pcolors),
	                     norm=cnorm) 

            #cm=ax.pcolormesh(axis[0,:], axis[1,:], pldat, cmap=my_cmap0, norm=cnorm) 
	else:
            cm=ax.pcolormesh(axis[0,:], axis[1,:], pldat, norm=cnorm) 


        if mask!=None:

            # get the color map of mask
	    """
            color1=colors.colorConverter.to_rgba('white')
            color2=colors.colorConverter.to_rgba('red')
            my_cmap=colors.LinearSegmentedColormap.from_list('mycmap',[color1, color2], 512) 
            my_cmap._init()
            alphas=np.linspace(0.2, 0.7, my_cmap.N+3)
            my_cmap._lut[:,-1] = alphas 
	    """
    
	    maskdata=np.ma.masked_where((mask<=1e-2)&(mask>=-1e-2) , mask)
            mymap=ax.contourf(axis[0,:], axis[1,:], maskdata, cmap=maskcolors)

            cbar=fig.colorbar(mymap, ticks=[4, 6, 8]) #, orientation='horizontal')
            cbar.ax.set_yticklabels(['void', 'filament', 'halo'])

	pl.xscale(axiscale[0])
	pl.yscale(axiscale[1])


        return 






    def hist1d(self, data, output, nbin=200, fact_xaxis=[1,1], label=['$log(1+\delta)$', 'PDF'], axiscale=['log', 'log']):
        ltype=mut.get_ltype()
        pl.clf()

        print 'hist1d: data shape: {0}'.format(data.shape)

	px=[]
        for i in range(data.shape[0]):
            py, edge=np.histogram(data[i], bins=nbin, density=True)
            px.append( np.delete(edge, nbin-1) )
            pl.plot(px[i], py, ltype[i], linewidth=2)
 
        min, max=np.min(px), np.max(px)
	pl.xlim(min*fact_xaxis[0], max*fact_xaxis[1])

	pl.xscale(axiscale[0])
	pl.yscale(axiscale[1])

        pl.xlabel(label[0])
        pl.ylabel(label[1])
        pl.savefig(output)
    









class myfig_fromfile():
    def __init__(self, fname):
	self.fname=fname

	"""         >>>get the 1D-data from list of fname[n]<<<
	    self.data: contains the list of 'n'data from fname[n], 
            each of which (self.data[i])[a][b] is then a 2D-array, where a
	    is the colum of the data, and b denotes the row of the data.
	"""
    def getdata(self):
	nfile = len(self.fname)
	self.data = [ list(zip(* np.genfromtxt(self.fname[i]) )) for i in range(nfile) ] 
	self.ncol = [ len(self.data[i]) for i in range(nfile)]


	""" >>> get the 1D plot  <<< """
    def plot(self, output, xlim, ylim, findex=0, allcol=True, colindex=[1], \
	    xscale='log', yscale='log'):
	pl.clf()
	if(allcol==True):
	    plots = range(self.ncol[findex]-1)
	else:
	    plots = list(np.array(colindex)-1)

	#for i in range(self.ncol[findex]-1):
	for i in plots:
	    pl.plot(self.data[findex][0], self.data[findex][i+1], 'k-', linewidth=1.3)

	pl.xscale(xscale)
	pl.yscale(yscale)
	pl.xlim(xlim)
	pl.ylim(ylim)

	pl.savefig(output)





""" >>>class 'mycontour' contains all methods to make 2 dimensional contour plots<<< """
class mycontour():
    def __init__(self, fname_data, fname_axis):
	self.fname_data=fname_data
        self.fname_axis=fname_axis 

    # get the 2D-data from list of fname[n]
    def getdata(self):
	nfile = len(self.fname_data)
	self.data, self.xaxis, self.yaxis = [], [], []
	for i in range(nfile):
	    dat=np.genfromtxt(self.fname_data[i])
	    data=np.array(list(zip(*dat)))
	    self.data.append(data)

            ax=np.genfromtxt(self.fname_axis[i]) 
	    xaxis=np.array(list(zip(*ax))[0])
	    yaxis=np.array(list(zip(*ax))[1])
	    self.xaxis.append(xaxis)
	    self.yaxis.append(yaxis)


    def cont2d(self, output, findex=0):
        pl.clf()
	nl_f, nl_l = 20, 20

	dat = self.data[findex] 
	dat = np.ma.masked_where(dat<=0, dat)

	
        contf =pl.contourf(self.xaxis[findex], self.yaxis[findex], dat, \
			  locator=ticker.LogLocator(), cmap=pl.cm.get_cmap('Blues'))

	""" different levels
        maxx, minn  = np.log10(np.max(dat) ), np.log10(1e-50)
        level = np.power(10,  np.arange(minn, maxx, (maxx-minn)/nl_f ) )
        contf =pl.contourf(self.xaxis[findex], self.yaxis[findex], dat, level, \
			 norm=colors.LogNorm() ,cmap=pl.cm.get_cmap('Blues'))
	# contour lines
        contl =pl.contour(self.xaxis[findex], self.yaxis[findex], dat, nl_l, \
		    	    colors='k', linewidth=3, alpha=0.5, hold = 'on')
	"""

	pl.xscale('log')
	pl.yscale('log')
	pl.savefig(output)









"""
#----------------------------------------------------------------#
fig=pl.figure()
r=2
#----------------------------------------------------------------#
ax1=fig.add_subplot(111)
ax1.set_xlim(3e-3, 0.5)
#ax1.set_ylim(0., 4)   
ax1.set_ylim(8e2, 2e5)   

#----------------------------------------------------------------#
dat=np.genfromtxt(fname_g1[0])
data=list(zip(*dat))
k_g1_b1 = np.array(data[0])
pl_g1_b1 = np.array(data[1])
pn_g1_b1 = np.array(data[4])

#--------------------------------#
dat=np.genfromtxt(fname_mc[0])
data=list(zip(*dat))
k_mc_b1 = np.array(data[0])
pl_mc_b1 = np.array(data[1])
pn_mc_b1nparraydata[4])



#----------------------------------------------------------------#
plin_g1_b1 = pl_g1_b1* (3.*(np.sin(k_g1_b1*2) - k_g1_b1*r*np.cos(k_g1_b1*r) ) /np.power(k_g1_b1*r,3) )
pn_mc_b1_full = np.interp(k_g1_b1, k_mc_b1, pn_mc_b1)


ax1.loglog(k_g1_b1, (pn_g1_b1+pn_mc_b1_full), 'k-', linewidth=1.8)
ax1.loglog(k_g1_b1, pn_g1_b1, 'k--', dashes=(15, 10), linewidth=1.8)
ax1.loglog(k_g1_b1, pn_mc_b1_full, 'k--', linewidth=1.8)


ax1.loglog(k_g1_b1, plin_g1_b1, 'k:', linewidth=1.8)
ax1.loglog(k_g1_b1, plin_g1_b1*1.5*1.5, 'k-.', linewidth=1.8)

#-------------------------------------------------------------------

ax1.set_xlabel('$k(Mpc/h)$', fontsize='large')
ax1.set_ylabel('$P(k)$', fontsize='large')
ax1.text(0.004, 1.2e5, '$z=0.5$', fontsize=18)

line1 = pl.Line2D(range(10), range(10), linestyle='-', color='k', linewidth=1.8)
line2 = pl.Line2D(range(10), range(10), linestyle='--', dashes=(15,10), color='k', linewidth=1.8)
line3 = pl.Line2D(range(10), range(10), linestyle='--', color='k', linewidth=1.8)
line4 = pl.Line2D(range(10), range(10), linestyle=':', color='k', linewidth=1.8)
line5 = pl.Line2D(range(10), range(10), linestyle='-.', color='k', linewidth=1.8)

l1=pl.legend((line1, line4,line5), ('$P_E(k)$', '$P_L(k)$', '$b_1^2 P_L(k)$'), \
             handlelength=3, bbox_to_anchor =  (0.96, 0.96) )

l2=pl.legend((line2,line3), ('$\Gamma_E^{(1)}$', '$\Gamma_E^{(2)}$'), \
             handlelength=3, bbox_to_anchor =  (0.3, 0.32) )

pl.gca().add_artist(l1)

#-------------------------------------------------------------------
#ax1.set_xlabel('$k(Mpc/h)$', fontsize='large')
#ax1.set_ylabel('$P_E(k)/( b^2_{eff} P_L(k))$', fontsize='large')
#ax1.text(0.04, 1.55, '$z=3$', fontsize=18)



pl.savefig('../result/rsd/fig/pk_rsd_2D.png')
pl.show()
"""

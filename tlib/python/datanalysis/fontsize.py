import pylab as M

def mySetFontSize(xlabel=0, ylabel=0, label=0, legend=0,latexfont=1):
    """
    set the ticklabel sizes
    """
    if label > 0:
        xlabel = label
        ylabel = label
        
    ax = M.gca()
    if (xlabel > 0) or (ylabel > 0):
        if xlabel > 0:
            lx = ax.get_xticklabels()
            M.setp(lx,fontsize = xlabel)
            M.setp(lx,fontname = 'Times New Roman')
        if ylabel > 0:
            ly = ax.get_yticklabels()
            M.setp(ly,fontsize = ylabel)
            M.setp(ly,fontname = 'Times New Roman')

    if legend > 0:
        leg = ax.get_legend()
        ltext = leg.get_texts()  # all the text.Text instance in the legend
        M.setp(ltext, fontsize=legend)    # the legend text fontsize
        if latexfont == 1:
            M.setp(ltext,fontname = 'Times')

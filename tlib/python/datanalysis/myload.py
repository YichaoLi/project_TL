import numpy as N

def myload(fileName,ncol=1,dtype='float',skipnums=0,sep=' '):
    inFile = open(fileName,'r')
    ret = N.fromfile(inFile,dtype=dtype,sep=sep)[skipnums:]
    inFile.close()
    nrow = ret.shape[0]
    nrow = nrow/ncol
    ret.resize((nrow,ncol))
    return ret

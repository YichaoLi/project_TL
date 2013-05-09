import numpy as np



def extrapolate(x, y, x_out, order=2, extra_type='linear'):
    ''' >> Extrapolate the data with polynomial upto `order` << '''

    #print x.shape, y.shape

    if extra_type=='linear':
        xx, yy, xx_out =x, y, x_out
    elif extra_type=='log':
        xx, yy, xx_out = np.log(x), np.log(y), np.log(x_out)
    else:
        raise Exception()


    p=np.polyfit(xx, yy, order)
    y_out=np.polyval(p, xx_out)

    if extra_type=='linear':
        return y_out
    elif extra_type=='log':
        return np.exp(y_out)

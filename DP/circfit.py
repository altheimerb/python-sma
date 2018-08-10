#based on https://scipy-cookbook.readthedocs.io/items/Least_Squares_Circle.html
#using optimize.leastsq
#update: use optiimize.least_squares so that i can add constraints TO DO
#given x,y data, find circle fit. 
#uses all the data; if trimming trace is desired, do that in the function call
import numpy
from numpy import sqrt
from scipy import optimize

def cfit(x,y,xcenguess, ycenguess, radguess, bounds=(-numpy.inf,numpy.inf)):
	###beta0=[xcenguess, ycenguess, radguess] #initial guess
	# print(kwargs)
	# if kwargs is not {}:
		# print('kwargs present')
		# for key, value in kwargs.iteritems():
			# if key == 'bounds':
				# bo = value
			# else:
				# bo = (-np.inf, np.inf) #no bounds specified
	# else:
		# bo = (-np.inf, np.inf) #no bounds specified
		# print("no keyword")
	#print bounds
	center_estimate = xcenguess, ycenguess
	#print(center_estimate)
	#print(x)
	#print(y)
	#center_2b, ier = optimize.leastsq(f_2b, center_estimate, args=(x,y), Dfun=Df_2b, col_deriv=True)
	#xc_2, yc_2 = center_2b
	center_2b  =optimize.least_squares(f_2b,center_estimate, Df_2b, args=(x,y), bounds=bounds)
	xc_2, yc_2 = center_2b.x
	#Ri_2       = calc_R(x,y,*center_2b)
	Ri_2       = calc_R(x,y,xc_2,yc_2)
	R_2        = Ri_2.mean()
	residu_2   = sum((Ri_2 - R_2)**2)/(x.size) #avg sq residual
	#print(x.size)
	#print(xc_2)
	#print(yc_2)
	#print(R_2)
	return [xc_2, yc_2, R_2, residu_2]
	
def calc_R(x, y , xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return sqrt((x-xc)**2 + (y-yc)**2)

def f_2b(c,x,y):
    """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(x,y,*c)
    return Ri - Ri.mean()
	
def Df_2b(c,x,y):
    """ Jacobian of f_2b
    The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
    xc, yc     = c
    df2b_dc    = numpy.empty((len(c), x.size))

    Ri = calc_R(x,y,xc, yc)
    df2b_dc[0] = (xc - x)/Ri                   # dR/dxc
    df2b_dc[1] = (yc - y)/Ri                   # dR/dyc
    df2b_dc    = df2b_dc - df2b_dc.mean(axis=1)[:, numpy.newaxis]

    return numpy.transpose(df2b_dc)
	
	
# def cirfn(beta,x)
	# """ implicit definition of the circle """
    # return (x[0]-beta[0])**2 + (x[1]-beta[1])**2 -beta[2]**2
	#beta[0] is xcenter, beta[1] y center, beta[2] radius
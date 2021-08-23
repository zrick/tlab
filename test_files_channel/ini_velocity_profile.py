import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
from scipy import integrate
# import matplotlib.colors as mcolors

#---------------------------------------------------------------------------#
# path to 3d-fields
path  = str(os.path.dirname(__file__) + '/../test_parallel_channel/' ) # name = 'flow.20.1' # file = str(path+name)
# index = 0

#---------------------------------------------------------------------------#
# read grid file
grid = mp.DnsGrid(path+'grid')

# # read flow field
# flow = mp.Field(path,var='flow',index=index)
# flow.read_3d_field()

#---------------------------------------------------------------------------#
# ProfileVelocity    = parabolic
VelocityX            = 0.0
YCoorVelocity        = 0.5#1.0	# y-coordinate of profile reference point, relative to the total scale, equation (2.1).
DeltaVelocity        = 1.0	# Reference profile difference, equation (2.1).

# compute ThickVelocity to ensure u(y=0)=0
ThickVelocity        = grid.y.max() * YCoorVelocity / (2 * np.sqrt(1 - VelocityX / DeltaVelocity))
# in this case: ThickVelocity = 0.5 * grid.y.max() * YCoorVelocity
print('ThickVelocity to ensure u(y=0)=0: ', str(ThickVelocity)) # Reference profile thickness, equation (2.1).

# compute u_profile
y_start = grid.y[0]
y_end   = grid.y[-1]
ycenter = y_start + y_end * YCoorVelocity
yrel    = grid.y - ycenter

xi = yrel / ThickVelocity
amplify = ( 1 - xi/2 ) * ( 1 + xi/2 )

u_profile = VelocityX + DeltaVelocity * amplify 

#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# plot u(y) -- ini profile
plt.figure(figsize=size)
plt.xlabel("u_mean-velocity")
plt.ylabel("y")
plt.xlim(0,1)
plt.ylim(0,grid.y.max())
plt.grid('True')
plt.plot(u_profile, grid.y, marker='.',label='u_mean')
plt.legend(loc=1)
plt.show()
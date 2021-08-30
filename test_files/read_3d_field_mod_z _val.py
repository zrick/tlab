# %%
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
from   scipy import interpolate
# import matplotlib.colors as mcolors

#---------------------------------------------------------------------------#
# path to 3d-fields
path  = str(os.path.dirname(__file__) + '/../test_parallel/' ) # name = 'flow.20.1' # file = str(path+name)
index = 0

#---------------------------------------------------------------------------#
# read grid and flow fields
grid = mp.DnsGrid(path+'grid')

# orignial field
flow = mp.Field(path,var='flow',index=index)
flow.read_3d_field()

# u_mod field 
f = open(path +'fld_mod.1','rb')
f.seek(52,0)
u_mod = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
u_mod = u_mod.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# read eps field
f = open(path +'eps0.1','rb')
f.seek(52,0)
eps = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
eps = eps.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# %%
#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
# 2d plot - yz
plt.figure(figsize=size)
plt.title('2d-plot -- yz-plane -- velocity u_mod')
plt.xlabel("z")
plt.ylabel("y")
#
plt.pcolormesh(grid.z,grid.y,u_mod[0,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xy
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u_mod')
plt.xlabel("x")
plt.ylabel("y")
#
plt.pcolormesh(grid.x,grid.y,u_mod[:,:,grid.nz//3].T, shading=shading, cmap='RdBu_r') #shading='gouraud',
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xz
plt.figure(figsize=size)
plt.title('2d-plot -- xz-plane -- velocity u_mod')
plt.xlabel("x")
plt.ylabel("z")
#
plt.pcolormesh(grid.x,grid.z,u_mod[:,1,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#

# %%

plt.figure(figsize=size)
plt.xlabel("z")
plt.ylabel("w-velocity")
for i in range(0,10):
    plt.plot(grid.z,u_mod[-1,i,:], marker='.',label='y-node='+str(i))
plt.legend(loc=1)
plt.show()

# small check 
w     = flow.w * (1. - eps)
w_mod = u_mod  * (1. - eps)
res   = w - w_mod
print(str(res.sum()))

# %%  validate splines
i = 5

plt.figure(figsize=size)
plt.xlabel("z")
plt.ylabel("w-velocity")
plt.plot(grid.z,u_mod[0,i,:], marker='.',label='y-node='+str(i))
plt.legend(loc=1)
plt.show()

yk1 = u_mod[0,i,8:12]
yk2 = u_mod[0,i,20:24]
yk  = np.concatenate((yk1,yk2))
xk1 = grid.z[8:12]
xk2 = grid.z[20:24]
xk  = np.concatenate((xk1,xk2))
x_new = grid.z[11:21]

# spline interpolation
tck = interpolate.splrep(xk, yk, s=0, k=3)

# new
yks = interpolate.splev(x_new, tck, der=0)

plt.figure(figsize=size)
plt.xlabel("z")
plt.ylabel("w-velocity")
plt.plot(grid.z[8:24],u_mod[0,i,8:24], marker='.',label='y-node='+str(i))
plt.plot(x_new,yks, marker='.',label='new')
plt.legend(loc=1)
plt.show()

res = yks - u_mod[0,i,11:21]
print(str(res.sum()))


# %%
# plt.figure(figsize=size)
# plt.xlabel("z")
# plt.ylabel("w-velocity")
# for i in range(0,5):
#     plt.plot(grid.z,u_mod_glob[-1,i,:], marker='.',label='y-node='+str(i))
#     plt.plot(grid.z,u_mod_loc[-1,i,:], marker='x',label='y-node='+str(i))
# plt.legend(loc=1)
# plt.show()

# %%
# plt.figure(figsize=size)
# plt.xlabel("z")
# plt.ylabel("w-velocity")
# for i in range(0,5):
#     plt.plot(grid.z,u_mod_glob[-1,i,:] - u_mod_loc[-1,i,:], marker='.',label='y-node='+str(i))
#     # plt.plot(grid.z,u_mod_loc[-1,i,:], marker='x',label='y-node='+str(i))
# plt.legend(loc=1)
# plt.show()

sys.exit()
#---------------------------------------------------------------------------#
# 2d plot - yz
plt.figure(figsize=size)
plt.title('2d-plot -- yz-plane -- velocity u')
plt.xlabel("z")
plt.ylabel("y")
#
plt.pcolormesh(grid.z,grid.y,flow.w[0,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xy
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("y")
#
plt.pcolormesh(grid.x,grid.y,flow.w[:,:,grid.nz//3].T, shading=shading, cmap='RdBu_r') #shading='gouraud',
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xz
plt.figure(figsize=size)
plt.title('2d-plot -- xz-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("z")
#
plt.pcolormesh(grid.x,grid.z,flow.w[:,1,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
plt.show()





# %%
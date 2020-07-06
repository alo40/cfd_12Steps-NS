"""
2D Non-Linear convection as given on: 
https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/08_Step_6.ipynb
Today's 13.06.2020, some day of the Quarentene
"""


# preamble
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# spatial and array parameters
# ------------------------------------------------------
Co = 0.5                            # Courant number, >0.5 leads to instability
nx = 101                            # elements of x
ny = 101                            # elements of y
Lx = 1                              # x-length
Ly = 1                              # y-length
dx = Lx / (nx-1)                    # x-spatial step
dy = Ly / (ny-1)                    # y-spatial step
x = np.linspace(-Lx/2, Lx/2, nx)    # nx x 1 array
y = np.linspace(-Ly/2, Ly/2, ny)    # ny x 1 array
u = np.ones((nx, ny))               # x-velocity as nx x ny array 
v = np.ones((nx, ny))               # y-velocity as nx x ny array


# initial conditions 
# ------------------------------------------------------
u[int(0.25 / dy):int(0.75 / dy+1), int(0.25 / dx):int(0.75 / dx+1)] = 8 # hat function
v[int(0.25 / dy):int(0.75 / dy+1), int(0.25 / dx):int(0.75 / dx+1)] = 2 # hat function


# figure(s) configuration
# ------------------------------------------------------
fig, axs = plt.subplots(1, 2) # horizontal array

cax = [None, None] # list init
cax[0] = axs[0].pcolormesh(x, y, u) 
cax[1] = axs[1].pcolormesh(x, y, v) 

txt = plt.gcf().text(.15, .90, '') # outside the axis


# time step acc. to the CLF stability criterion
# ------------------------------------------------------
def timeStepCalc(cx, cy):
    
    # x-time step
    if cx == 0:
        dtx = 1000 # very big, so it is not selected
    else:
        dtx = Co * dx / cx

    # y-time step
    if cy == 0:
        dty = 1000 # very big, so it is not selected
    else:
        dty = Co * dy / cy      

    return min(dtx, dty)      # total time step, choose min value to assure numerical stability


# 2D stencil using array operations
# ------------------------------------------------------
def stencil2D(u, v):
    
    # copy u
    un = u.copy()
    vn = v.copy()
    
    # time step
    cx = np.max(u) # x max propagation speed
    cy = np.max(v) # y max propagation speed
    dt = timeStepCalc(cx, cy)
    
    # 2D Array operations
    # u(j  , i  ) --> u[1:-1, 1:-1]
    # u(j  , i-1) --> u[1:-1, 0:-2]
    # u(j-1, i  ) --> u[0:-2, 1:-1]
    
    # u stencil
    u[1:-1, 1:-1] = (un[1:-1, 1:-1] - (un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, 0:-2])) -
                                      (vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[0:-2, 1:-1])))

    # v stencil
    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - (un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2])) -
                                      (vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1])))

    # u boundary condtion (cyclic)
    u[ 0,  :] = u[-1,  :]
    u[-1,  :] = u[-2,  :]
    u[ :,  0] = u[ :, -1]
    u[ :, -1] = u[ :, -2]
    
    # v boundary condtion (cyclic)
    v[ 0,  :] = v[-1,  :]
    v[-1,  :] = v[-2,  :]
    v[ :,  0] = v[ :, -1]
    v[ :, -1] = v[ :, -2]
    
    return u, v # simple as that!


# main function
# ------------------------------------------------------
def animate(n):
    
    # update main text 
    txt.set_text('step = %4d' % (n))
    
    # declare u as global to be able to modify it
    global u, v
    
    # use 2D stencil
    u, v = stencil2D(u, v)

    # module of velocity vector (u, v)          
    U = np.sqrt(u**2 + v**2)    
    
    # update samples
    cax[0].set_array(u[:-1, :-1].flatten()) # [:-1, :-1] is to correct visualization bug
    cax[1].set_array(U[:-1, :-1].flatten()) # [:-1, :-1] is to correct visualization bug
    
    
# animation
# ------------------------------------------------------
anim = FuncAnimation(fig, animate, frames=100, interval=100, repeat=False)
# anim.save('./NavierStokes_Step05.gif', writer='imagemagick', fps=30, dpi=100) # dpi>100 is to cumbersome


# show maximized figure
# ------------------------------------------------------
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized() 
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()


# not used
# ------------------------------------------------------
# print("n = %d, m = %d" % (n, m)) # for debugging only

"""
Poisson equation as given on: 
https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/13_Step_10.ipynb
Today's 20.06.2020, this one's like Laplace but with a source term!
"""

# preamble
# ------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.cm as cm
import myFunctions as mi


# user defined parameters
# ------------------------------------------------------
Co = 0.1                            # Courant number, > 0.5 leads to instability in convection, > 0.1 in diffusion
Nu = 0.1                            # diffusion coefficient
nx = 81                             # elements of x
ny = 81                             # elements of y
Lx = 1                              # x-axis length
Ly = 1                              # y-axis length
hmin = 0.35                         # IC hat square min size
hmax = 0.65                         # IC hat square max size
U = 10                              # IC hat u-velocity
V = 10                              # IC hat v-velocity
B = 1000                             # source term +/- max value


# arrays generation
# ------------------------------------------------------
dx = Lx / (nx-1)                    # x-spatial step
dy = Ly / (ny-1)                    # y-spatial step
dt = 0                              # time step init
x = np.linspace(0, Lx, nx)          # nx x 1 array
y = np.linspace(0, Ly, ny)          # ny x 1 array
u = np.ones((nx, ny))*1             # x-velocity as nx x ny array 
v = np.ones((nx, ny))*1             # y-velocity as nx x ny array
p = np.ones((nx, ny))*0             # used for Laplace stencil (pressure?)
b = np.ones((nx, ny))*0             # soruce term


# initial conditions (also set max min colorbar values)
# ------------------------------------------------------
u[int(hmin / dy): int(hmax / dy + 1), int(hmin / dx): int(hmax / dx + 1)] = U # hat function
v[int(hmin / dy): int(hmax / dy + 1), int(hmin / dx): int(hmax / dx + 1)] = V # hat function


# source term init
# ------------------------------------------------------
b[int(ny / 4), int(nx / 4)] = B 
b[int(ny * 3 / 4), int(nx * 3 / 4)] = -B


# figure(s) configuration
# ------------------------------------------------------
cax = [None, None] # list init

# many axis
fig, axs = plt.subplots(1, 2) 

cax[0] = axs[0].pcolormesh(x, y, p, cmap=cm.RdBu) 
cax[1] = axs[1].pcolormesh(x, y, p, cmap=cm.PiYG) 

fig.colorbar(cax[0], ax=axs[0], orientation="horizontal", fraction=0.025, extend='both')
fig.colorbar(cax[1], ax=axs[1], orientation="horizontal", fraction=0.025, extend='both')

txt = plt.gcf().text(.15, .90, '') # outside the axis


# main function
# ------------------------------------------------------
def animate(n):
        
    # declare u as global to be able to modify it
    global u, v, p
    
    # # max propagation speed
    # cx = np.max(u)  
    # cy = np.max(v) 
    
    # # 2D Burgers' equation
    # dt = mi.timeStep_Burgers(cx, cy, dx, dy, Co, Nu)
    # u, v = mi.stencil_Burgers(u, v, dx, dy, dt, Nu)

    # # 2D steady-state equations
    # p, L = mi.stencil_Laplace(p, dx, dy) # Laplace
    p, L = mi.stencil_Poisson(p, b, dx, dy) # Poisson, with source term
    if L < 1e-8 and n > 1:
        return # stop function

    # boundary conditions
    # u = mi.boundaries_constant(u, 0)
    # v = mi.boundaries_constant(v, 0)
    p = mi.boundaries_constant(p, 0)
    # u = mi.boundaries_derivative(u)
    # u = mi.boundaries_cyclic(u) % not very realistic
    # p = mi.boundaries_Laplace(p, y)
    
    # update samples    
    cax[0].set_array(p[:-1, :-1].flatten()) 
    cax[1].set_array(p[:-1, :-1].flatten()) 
    
    # update main text 
    # txt.set_text('step = %4d' % (n))
    # txt.set_text('n = %4d \ndt = %.8f' % (n, dt))
    txt.set_text('n = %4d \nL = %.8f' % (n, L))
    
    
# animation
# ------------------------------------------------------
anim = FuncAnimation(fig, animate, frames=1000, interval=10, repeat=False)


# show maximized figure
# ------------------------------------------------------
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized() 
plt.show()


# # 3D projection (only final step)
# # ------------------------------------------------------
# from mpl_toolkits.mplot3d import Axes3D
# X, Y = np.meshgrid(x, y)
# fig3d = plt.figure(figsize=(11, 7), dpi=100)
# axs3d = fig3d.gca(projection='3d')
# surf = axs3d.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.viridis, linewidth=0, antialiased=False)


# # for array testing
# # ------------------------------------------------------
# a = np.transpose(np.array([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]))
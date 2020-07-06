"""
Laplace equation as given on: 
https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/12_Step_9.ipynb
Today's 18.06.2020, quarentene may be over soon, but also this course!
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
nx = 41                             # elements of x
ny = 41                             # elements of y
Lx = 1                              # x-axis length
Ly = 1                              # y-axis length
hmin = 0.35                         # IC hat square min size
hmax = 0.65                         # IC hat square max size
U = 10                              # IC hat u-velocity
V = 10                              # IC hat v-velocity


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


# initial conditions (also set max min colorbar values)
# ------------------------------------------------------
u[int(hmin / dy): int(hmax / dy + 1), int(hmin / dx): int(hmax / dx + 1)] = U # hat function
v[int(hmin / dy): int(hmax / dy + 1), int(hmin / dx): int(hmax / dx + 1)] = V # hat function
# p[0, 0] = 10 # only as reference for the colorbar max value


# figure(s) configuration
# ------------------------------------------------------
cax = [None, None] # list init

# # one axis
# fig, axs = plt.subplots()
# cax[0] = axs.pcolormesh(x, y, u, cmap=cm.hot)
# fig.colorbar(cax[0], ax=axs)

# many axis
fig, axs = plt.subplots(1, 2) 

cax[0] = axs[0].pcolormesh(x, y, u, cmap=cm.hot_r) 
cax[1] = axs[1].pcolormesh(x, y, v, cmap=cm.gray_r) 

fig.colorbar(cax[0], ax=axs[0], orientation="horizontal", fraction=0.025)
fig.colorbar(cax[1], ax=axs[1], orientation="horizontal", fraction=0.025)

txt = plt.gcf().text(.15, .90, '') # outside the axis


# main function
# ------------------------------------------------------
def animate(n):
    
    # declare u as global to be able to modify it
    global u, v, p
    
    # max propagation speed
    cx = np.max(u)  
    cy = np.max(v) 
    
    # # 2D convection
    # dt = mi.timeStep_convection(cx, cy, dx, dy, Co)
    # u, v = mi.stencil_convection(u, v, dx, dy, dt)
    
    # # 2D diffusion (only u)
    # dt = mi.timeStep_diffusion(dx, dy, Co, Nu)
    # u = mi.stencil_diffusion(u, dx, dy, dt, Nu)

    # # 2D Burgers' equation
    # dt = mi.timeStep_Burgers(cx, cy, dx, dy, Co, Nu)
    # u, v = mi.stencil_Burgers(u, v, dx, dy, dt, Nu)

    # # 2D Laplace equation
    # p, L = mi.stencil_Laplace(p, dx, dy)
    # if L < 1e-8 and n > 1:
    #     return # stop function

    # boundary conditions
    # u, v = mi.boundaries_constant(u, v, 0)
    u, v = mi.boundaries_derivative(u, v)
    # u, v = mi.boundaries_cyclic(u, v) % not very realistic
    # p = mi.boundaries_Laplace(p, y)
    
    # update samples
    # cax.set_array(p[:-1, :-1].flatten()) # [:-1, :-1] is to correct visualization bug
    
    cax[0].set_array(u[:-1, :-1].flatten()) 
    cax[1].set_array(v[:-1, :-1].flatten()) 
    
    # update main text 
    # txt.set_text('step = %4d' % (n))
    txt.set_text('n = %4d \ndt = %.8f' % (n, dt))
    # txt.set_text('n = %4d \nL = %.8f' % (n, L))
    
    
# animation
# ------------------------------------------------------
anim = FuncAnimation(fig, animate, frames=1000, interval=25, repeat=False)
# anim.save('./NavierStokes_Step07.gif', writer='imagemagick', fps=30, dpi=100) # dpi>100 is to cumbersome


# show maximized figure
# ------------------------------------------------------
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized() 
# plt.gca().get_visibleset_aspect('equal', adjustable='box')
plt.show()


# for array testing
# ------------------------------------------------------
a = np.transpose(np.array([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]))
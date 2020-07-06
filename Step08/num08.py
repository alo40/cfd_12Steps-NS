"""
Burgers' equation as given on: 
https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/10_Step_8.ipynb
Today's 16.06.2020, I probably should be sleeping ...
"""

# preamble
# ------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import myFunctions as mi


# user defined parameters
# ------------------------------------------------------
Co = 0.1                            # Courant number, > 0.5 leads to instability in convection, > 0.1 in diffusion
Nu = 0.1                            # diffusion coefficient
nx = 101                            # elements of x
ny = 101                            # elements of y
Lx = 1                              # x-length
Ly = 1                              # y-length


# arrays generation
# ------------------------------------------------------
dx = Lx / (nx-1)                    # x-spatial step
dy = Ly / (ny-1)                    # y-spatial step
dt = 0                              # time step init
x = np.linspace(-Lx/2, Lx/2, nx)    # nx x 1 array
y = np.linspace(-Ly/2, Ly/2, ny)    # ny x 1 array
u = np.ones((nx, ny))*0             # x-velocity as nx x ny array 
v = np.ones((nx, ny))*0             # y-velocity as nx x ny array


# initial conditions 
# ------------------------------------------------------
u[int(0.40 / dy):int(0.60 / dy+1), int(0.40 / dx):int(0.60 / dx+1)] = 24 # hat function
v[int(0.40 / dy):int(0.60 / dy+1), int(0.40 / dx):int(0.60 / dx+1)] = 24 # hat function


# figure(s) configuration
# ------------------------------------------------------
fig, axs = plt.subplots(1, 2) # horizontal array

cax = [None, None] # list init
cax[0] = axs[0].pcolormesh(x, y, u) 
# cax[1] = axs[1].pcolormesh(x, y, v) 

txt = plt.gcf().text(.15, .90, '') # outside the axis


# main function
# ------------------------------------------------------
def animate(n):
    
    # declare u as global to be able to modify it
    global u, v
    
    # time step
    cx = np.max(u) # x max propagation speed
    cy = np.max(v) # y max propagation speed
    
    # # 2D convection
    # dt = mi.timeStep_convection(cx, cy, dx, dy, Co)
    # u, v = mi.stencil_convection(u, v, dx, dy, dt)
    
    # # 2D diffusion (only u)
    # dt = mi.timeStep_diffusion(dx, dy, Co, Nu)
    # u = mi.stencil_diffusion(u, dx, dy, dt, Nu)

    # 2D Burgers' equation
    dt = mi.timeStep_Burgers(cx, cy, dx, dy, Co, Nu)
    u, v = mi.stencil_Burgers(u, v, dx, dy, dt, Nu)

    # boundary conditions
    u, v = mi.boundaries_cyclic(u, v)

    # # module of velocity vector (u, v)          
    # U = np.sqrt(u**2 + v**2)    
    
    # update samples
    cax[0].set_array(u[:-1, :-1].flatten()) # [:-1, :-1] is to correct visualization bug
    # cax[1].set_array(U[:-1, :-1].flatten()) # [:-1, :-1] is to correct visualization bug
    
    # update main text 
    # txt.set_text('step = %4d' % (n))
    txt.set_text('n = %4d \ndt = %.8f' % (n, dt))
    
    
# animation
# ------------------------------------------------------
anim = FuncAnimation(fig, animate, frames=500, interval=10, repeat=False)
# anim.save('./NavierStokes_Step07.gif', writer='imagemagick', fps=30, dpi=100) # dpi>100 is to cumbersome


# show maximized figure
# ------------------------------------------------------
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized() 
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()

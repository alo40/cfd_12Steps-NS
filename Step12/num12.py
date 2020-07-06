"""
Poisson equation as given on: 
https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/15_Step_12.ipynb
Today's 05.07.2020, last step!! This ends today .... 
"""

# preamble
# ------------------------------------------------------
from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib as mpl
import matplotlib.cm as cm
import myFunctions as mi


# user defined parameters
# ------------------------------------------------------
Co = 0.1                            # Courant number, > 0.5 leads to instability in convection, > 0.1 in diffusion
Nu = 0.1                            # diffusion coefficient
rho = 1                             # density
nx = 401                            # elements of x
ny = 401                            # elements of y
Lx = 1                              # x-axis length
Ly = 1                              # y-axis length
hmin = 0.4                          # IC hat square min size
hmax = 0.6                          # IC hat square max size
U = 0                               # IC hat u-velocity
V = 0                               # IC hat v-velocity
B = 1                               # source term +/- max value


# arrays generation
# ------------------------------------------------------
dx = Lx / (nx-1)                    # x-spatial step
dy = Ly / (ny-1)                    # y-spatial step
dt = 0                              # time step init
x = np.linspace(0, Lx, nx)          # nx x 1 array
y = np.linspace(0, Ly, ny)          # ny x 1 array
u = np.ones((nx, ny))*0             # x-velocity as nx x ny array 
v = np.ones((nx, ny))*0             # y-velocity as nx x ny array
p = np.ones((nx, ny))*0             # used for Laplace stencil (pressure?)
b = np.ones((nx, ny))*B             # source term
X, Y = np.meshgrid(x, y)            # only for the vector field
m = int(nx/40)                      # only for the vector field, show every m values


# initial conditions (also set max min colorbar values)
# ------------------------------------------------------
# u[0, 0] = 1
# v[0, 0] = 1
# p[0, 0] = 10
# p[1, 0] = -10
# u[int(hmin / dy): int(hmax / dy + 1), int(hmin / dx): int(hmax / dx + 1)] = U # hat function
# v[int(hmin / dy): int(hmax / dy + 1), int(hmin / dx): int(hmax / dx + 1)] = V # hat function


# source term init
# ------------------------------------------------------
# b[int(ny / 4), int(nx / 4)] = B 
# b[int(ny * 3 / 4), int(nx * 3 / 4)] = -B


# figure(s) configuration
# ------------------------------------------------------
cax = [None, None]
fig, axs = plt.subplots(1, 2)

# colormesh plot
cax[0] = axs[0].pcolormesh(x, y, u, cmap=cm.hot)
# cax[1] = axs[1].pcolormesh(x, y, v, cmap=cm.hot)
cax[1] = axs[1].pcolormesh(x, y, p, cmap=cm.PiYG)

# quiver plot
# cax[1] = axs[1].quiver(X[::m, ::m], Y[::m, ::m], u[::m, ::m], v[::m, ::m], pivot='mid', color='r', 
#                         scale=1/0.1,
#                         headwidth=0,
#                         headlength=0,
#                         )

fig.colorbar(cax[0], ax=axs[0], orientation="horizontal", fraction=0.025)
fig.colorbar(cax[1], ax=axs[1], orientation="horizontal", fraction=0.025)#, extend='both')

txt = plt.gcf().text(.15, .90, '') # outside the axis


# main function
# ------------------------------------------------------
def animate(n):
        
    # declare u as global to be able to modify it
    global u, v, p
    
    # max propagation speed
    cx = np.max(u)  
    cy = np.max(v) 
       
    # 2D Navier-Stokes equations
    dt = mi.timeStep_NavierStokes(cx, cy, dx, dy, Co, Nu, rho)
    u, v, p = mi.stencil_NavierStokes(u, v, p, b, dx, dy, dt, Nu, rho)

    # # constant boundary conditions (e.g. u = a)
    u = mi.boundaries_constant(u, 0)
    v = mi.boundaries_constant(v, 0)
    # p = mi.boundaries_constant(p, 0)
    
    # derivative boundary conditions (e.g. du/dx = 0)
    # u = mi.boundaries_derivative(u)
    # v = mi.boundaries_derivative(v)
    p = mi.boundaries_derivative(p)
    
    # extra boundary conditions
    p[-1, :] = 0    # p = 0 at y = 2
    u[-1, :] = 1    # set velocity on cavity lid equal to 1
    
    # vector norm 
    u_norm = np.sqrt(u**2 + v**2) 
    
    # update samples  
    cax[0].set_array(u_norm[:-1, :-1].flatten()) # colormesh plot update
    cax[1].set_array(p[:-1, :-1].flatten()) # colormesh plot update
    # cax[1].set_UVC(u[::m, ::m], v[::m, ::m]) # quiver plot update
    
    # steamline plot
    # if n > 0 and n % 1000 == 0:
        # axs[1].collections = [] # clear lines streamplot
        # axs[1].patches = [] # clear arrowheads streamplot
        # cax[1] = axs[1].streamplot(X, Y, u, v, density=4, cmap='jet',arrowsize=1) # streamline plot update
    
    # update main text 
    txt.set_text('n = %4d \ndt = %.8f' % (n, dt))
    
    
# animation
# ------------------------------------------------------
anim = FuncAnimation(fig, animate, frames=1001, interval=250, repeat=False)


# show maximized figure
# ------------------------------------------------------
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized() 
plt.show()


# 2D quiver plot (only final step)
# ------------------------------------------------------
"""
m = 16
fig, ax = plt.subplots(1,1)
X, Y = np.meshgrid(x, y)
ax.quiver(X[::m,::m], Y[::m,::m], u[::m,::m], v[::m,::m], pivot='mid', color='r', units='inches')
plt.streamplot(X, Y, u, v)
"""


# 3D projection (only final step)
# ------------------------------------------------------
"""
from mpl_toolkits.mplot3d import Axes3D
X, Y = np.meshgrid(x, y)
fig3d = plt.figure(figsize=(11, 7), dpi=100)
axs3d = fig3d.gca(projection='3d')
surf = axs3d.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.viridis, linewidth=0, antialiased=False)
"""


# for array testing (not used)
# ------------------------------------------------------
a = np.transpose(np.array([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]))


# 2x2 subplots (not used, too small)
# ------------------------------------------------------
# cax = [None, None, None, None] # list init
# fig, axs = plt.subplots(2, 2) 

# cax[0] = axs[0, 0].pcolormesh(x, y, u, cmap=cm.hot) 
# cax[1] = axs[0, 1].pcolormesh(x, y, v, cmap=cm.hot)
# cax[2] = axs[1, 0].pcolormesh(x, y, p, cmap=cm.PiYG) 
# cax[3] = axs[1, 1].pcolormesh(x, y, p, cmap=cm.PiYG)

# fig.colorbar(cax[0], ax=axs[0, 0], orientation="vertical", fraction=0.025)
# fig.colorbar(cax[1], ax=axs[0, 1], orientation="vertical", fraction=0.025)
# fig.colorbar(cax[2], ax=axs[1, 0], orientation="vertical", fraction=0.025, extend='both')
# fig.colorbar(cax[3], ax=axs[1, 1], orientation="vertical", fraction=0.025, extend='both')


# other stencils (not used)
# ------------------------------------------------------
   # # 2D Burgers' equation
    # dt = mi.timeStep_Burgers(cx, cy, dx, dy, Co, Nu)
    # u, v = mi.stencil_Burgers(u, v, dx, dy, dt, Nu)

    # # 2D steady-state equations
    # p, L = mi.stencil_Laplace(p, dx, dy) # Laplace
    # p, L = mi.stencil_Poisson(p, b, dx, dy) # Poisson, with source term
    # if L < 1e-8 and n > 1:
        # return # stop function


# streamline plot
# ------------------------------------------------------
# config
# cax[1] = axs[1].streamplot(X, Y, u, v)
# axs[1].set_xlim([0, 1])
# axs[1].set_ylim([0, 1]) 

    # update
    # axs[1].collections = [] # clear lines streamplot
    # axs[1].patches = [] # clear arrowheads streamplot
    # cax[1] = axs[1].streamplot(X, Y, u, v, density=4, cmap='jet',arrowsize=1) # streamline plot update
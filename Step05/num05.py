"""
2D-Linear convection as given on: 
https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/07_Step_5.ipynb 
Today's 25.05.2020, some day of the Quarentene
"""

# preamble
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# spatial var
nx = 101                # elements of x
ny = 101                # elements of y
Lx = 2                  # x-length
Ly = 2                  # y-length

# numerical var
dx = Lx / (nx-1)        # x-spatial step
dy = Ly / (ny-1)        # y-spatial step
cx = 1.0                # x-propagation speed
cy = 1.0                # y-propagation speed    
Co = 0.1                # Courant number

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

dt = min(dtx, dty)      # total time step, choose min value to assure numerical stability

# arrays
x = np.linspace(-Lx/2, Lx/2, nx)    # nx x 1
y = np.linspace(-Ly/2, Ly/2, ny)    # ny x 1
u = np.zeros((nx, ny))              # nx x ny

# IC (hat!)
# u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2 
u[int(0.0 / dy):int(0.5 / dy + 1),int(0.0 / dx):int(0.5 / dx + 1)] = 2

# figure
fig, ax = plt.subplots()
cax = ax.pcolormesh(x, y, u)
fig.colorbar(cax)
txt = plt.gcf().text(.15, .90, '') # outside the axis

# main function
def animate(n):
    
    un = u.copy()   # update un
    
    # Array operations
    # u(i+1) --> u[2:]
    # u(i)   --> u[1:-1]
    # u(i-1) --> u[:-2]
    u[1:, 1:] = (un[1:, 1:] - (cx * dt / dx * (un[1:, 1:] - un[1:, :-1])) -
                              (cy * dt / dy * (un[1:, 1:] - un[:-1, 1:])))
    
    # BC
    u[0, :] = 0
    u[-1, :] = 0
    u[:, 0] = 0
    u[:, -1] = 0
    
    # update values
    txt.set_text('step = %4d \ndt = %.8f' % (n, dt))
    # frame_text.set_text('frame = %s \ntime = %.6s \ndt = %.6s' % (n, t, dt))
    cax.set_array(u[:-1, :-1].flatten()) # U given as an 1D array, the [:-1, :-1] is to correct visualization bug

# animation
anim = FuncAnimation(fig, animate, frames=1000, interval=100, repeat=False)
# anim.save('./NavierStokes_Step05.gif', writer='imagemagick', fps=30, dpi=100) # dpi>100 is to cumbersome
plt.draw()
plt.show()


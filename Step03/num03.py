# I'll try to make the code from scratch, taking inspiration from the given samples of the numpy page
# Today is th 8th of April 2020 (day 24th of the Quarantine)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

## space discretization
Lx = 8 # x-domain length 
nx = Lx*101 # number of elements 
dx = Lx/(nx-1) # space step

## time discretization & others
Co = 0.5 # Courant number, should be =< 0.5 for numerical stability
v = 0.1 # difussion coefficient
dt = Co*dx**2/v # time step
T = 1 # total time (sec?)
nt = round(T/dt) # number of frames

class UpdateDist(object):

    ## class init function
    def __init__(self, ax):
        self.line, = ax.plot([], [], lw=2)
        self.ax = ax
        self.ax.set_xlim(-Lx/2, Lx/2)
        self.ax.set_ylim(-0.2, 1.2)
        self.ax.grid(True)
        
        ## text parameters
        self.frame_text = ax.text(0.05, 0.1, '', transform=ax.transAxes)
        self.time_text = ax.text(0.05, 0.05, '', transform=ax.transAxes)
        
        ## values definition
        self.x = np.linspace(-Lx/2, Lx/2, nx)
        conds = [(self.x >= -Lx/2)  & (self.x <= -Lx/12),
                 (self.x >= -Lx/12) & (self.x <=  Lx/12),
                 (self.x >=  Lx/12) & (self.x <=  Lx/2)]
        funcs = [ 0., 0.5, 0.]
        self.u = np.piecewise(self.x, conds, funcs)
        self.u0 = self.u.copy()
       
    ## animation init function
    def init(self):
        self.line.set_data([], [])
        return self.line,

    ## main animation loop
    def __call__(self, n):

        ## require to continuously call the function if repeat=True
        if n == 0:
            self.u = self.u0.copy()
            return self.init()

        ## numerical scheme 
        un = self.u.copy()

        # ## using a for-loop
        # for i in range(len(self.x)-1):
        #     self.u[i] = v*dt/dx**2*(un[i+1] - 2*un[i] + un[i-1]) + un[i]

        ## using array operations (much faster)
        ## u(i+1) --> u[2:]
        ## u(i)   --> u[1:-1]
        ## u(i-1) --> u[:-2]
        self.u[1:-1] = v*dt/dx**2*(un[2:] - 2*un[1:-1] + un[:-2]) + un[1:-1]
            
        ## boundary conditions
        self.u[0]    = 0.
        self.u[nx-1] = 0.
        # self.u[0] = self.u[1]
        # self.u[nx-1] = self.u[nx-2]
            
        ## update info and return lines
        self.frame_text.set_text('frame = %4d' % (n))
        self.time_text.set_text('time = %2.4f' % (n*dt))
        self.line.set_data(self.x, self.u)
        return self.line, self.frame_text, self.time_text,

fig, ax = plt.subplots()
ud = UpdateDist(ax)
anim = FuncAnimation(fig, ud, frames=nt, init_func=ud.init,
                     repeat=False, interval=1, blit=True)
plt.show()

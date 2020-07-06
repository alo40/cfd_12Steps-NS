## Step 04 Burgers' Equation
## I know reusing code is a nice thing, but I need to learn and I only learn by doing
## Start on the 11th April 2020 (27th of the Quarantine)

## preamble
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

## space discretization
Lx = 4. # length of the x-domain, declared as float
nx = 101 # number of elements of th x-domain
dx = Lx/(nx-1) # space step
x = np.linspace(-Lx/2, Lx/2, nx) # x as numpy array

## time discretization
Co1 = 100.0 # Courant number for non-linear convection
Co2 = 0.1 # Courant number for diffusion
dt = 0. # init time step as float
c = 2.0 # constant propagation speed (only for testing)
v = 0.0001 # diffusion speed

## init function (hat!)
u = np.zeros(x.size)
conds = [(x >= -Lx/2) & (x <= -Lx/4),
         (x >= -Lx/4) & (x <=  Lx/4),
         (x >=  Lx/4) & (x <=  Lx/2)]
funcs = [0., 0.8, 0.]
u = np.piecewise(x, conds, funcs)

## figure setup
fig, ax = plt.subplots() # declare both fig and ax in one command
ax.set_xlim(-Lx/2, Lx/2) # x-axis limits acc. to Lx
# ax.set_ylim(np.min(u) - 0.1, np.max(u) + 0.1) # y-axis limits acc. to min/max value of u
ax.set_ylim(-0.1, 1.0)
plt.grid()

## figure outputs
line, = ax.plot([], [], lw=2) # if comma is missing, line is declared as 'list' and not as matplotlib object
text = ax.text(0.025, 0.75, '', transform=ax.transAxes) # information text

## animation!
def animate(n):

    ## init
    Convection = 0.0
    Diffusion = 0.0
    un = u.copy()

    ## time step calculation
    dt1 = dx*Co1/np.max(np.absolute(u)) # for convection scheme
    dt2 = dx**2*Co2/v # for diffusion scheme
    dt = min(dt1, dt2) # choose min value to assure numerical stability 
    
    ## array operators!
    ## u(i+1) --> u[2:]
    ## u(i)   --> u[1:-1]
    ## u(i-1) --> u[:-2]
   
    # ## non-linear convection scheme
    # Convection = - un[1:-1]*dt/dx*(un[1:-1] - un[:-2])
    
    # ## diffusion scheme
    Diffusion = v*dt/dx**2*(un[2:] - 2*un[1:-1] + un[:-2])
    
    ## main scheme
    u[1:-1] = un[1:-1] + Convection + Diffusion
    
    # ## boundary conditions (fix)
    # u[0] = 0.
    # u[-1] = 0.
    
    ## boundary conditions (du/dx = 0)
    u[0] = u[1]
    u[-1] = u[-2]
    
    line.set_data(x, u)
    text.set_text('step = %4d' % (n))
    # text.set_text('n = %4d \ndt = %.8f \ndt1 = %.8f \ndt2 = %.8f' % (n, dt, dt1, dt2)) % only for testing
    return line, text,

anim = FuncAnimation(fig, animate, frames=2500, interval=1, repeat=False, blit=True)
# anim.save('NavierStokes_Step04.gif', writer='imagemagick', fps=30, dpi=75) # dpi>100 is to cumbersome
plt.show() # comment anim.save to see the figure properly

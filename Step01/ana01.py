# Anaytical solution to the 'Step 1: 1-D Linear Convection' problem
# Animated, done  on the 7th August 2019

## preamble
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

## geometry setup (same as num01.py)
Lx = 25.0 # length given as float
nx = 1001 # 1D mesh size
dx = 2*Lx/(nx-1) # space step size
Ly = 2 # wave max amplitude + a bit more

## figure setup (same as num01.py)
fig = plt.figure()
ax = plt.axes(xlim=(-Lx, Lx), ylim=(-1, Ly + 0.5)) # working
line, = ax.plot([], [], lw=2)
frame_text = ax.text(0.01, 0.95, '', transform=ax.transAxes) # for printing text in the plot

## wave setup
# c = np.ones(nx)*0.1 # wave velocity, why is it an array?
c = 1
x = np.linspace(-Lx, Lx, nx) # 1D x-space
u0 = np.zeros(nx) # init in 1
# u0[int(0.5/dx):int(1.5/dx)] = 2 # hat functionxs
# u0 = np.exp(-x) # exponential function
u0 = np.sin(2*np.pi*x/Lx) # sin function, period = Lx
# u0 = 1/(1 + np.exp(-x)) # sigmoid function
u = u0.copy() 
plt.ylim(np.min(u)-0.1, np.max(u)+0.1) # resize of the y-axis limits 

## initialization functio, use only to restart values, not actually displayed 
def init():
    line.set_data([], []) # empty set
    frame_text.set_text('') # empty text
    return line, frame_text 

def animate(n): # for n in range(nt)
    if n != 0: # first value equal to u0
        for m in range(len(x)):
            u[m] = u0[m - c*n] 
    
    line.set_data(x, u)
    frame_text.set_text('frame = %s' % (n))
    return line, frame_text

anim = animation.FuncAnimation(fig, animate, frames=1000, interval=10, repeat = False,
                               init_func=init, blit=True)
    
## plot setup and run!
plt.title('Step 1: 1-D Linear Convection (Analytical Solution)')
plt.xlabel('x-axis')
plt.ylabel('wave amplitude')
plt.show()

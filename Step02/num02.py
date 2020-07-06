## Step 02 Non-linear convection
## From: https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/02_Step_2.ipynb
## Start on the 8th September 2019
## Update on the 28th March 2020 (Quarantine!)

## preamble
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

## figure setup
fig = plt.figure()
ax = plt.axes() # no limits defined
line, = ax.plot([], [], lw=2)
frame_text = ax.text(0.5, 0.9, '', transform=ax.transAxes) # for printing text in the plot

## space discretization
Lx = 2. # length as float
nx = 1001 # mesh elements
dx = Lx/(nx-1) # space step, needed for the CLF condition1
x = np.linspace(-Lx, Lx, nx) # space as numpy array

## time discretization
dt = 0. # time step init as float
Co = 1. # CLF stability condition
# time_text = ax.text(0.5, 0.9, '', transform=ax.transAxes) # for printing time in the plot

## picewise function using a numpy method!
u = np.zeros(nx) # init as zeros
## Stairs function!
# conds = [(x >= -2) & (x <= -1),
#          (x >= -1) & (x <=  0),
#          (x >=  0) & (x <=  1),
#          (x >=  1) & (x <=  2)]
# funcs = [ 0.1, 0.2, 0.3, 0.4]
## Hat function!
conds = [(x >= -2) & (x <= -1),
         (x >= -1) & (x <=  1),
         (x >=  1) & (x <=  2)]
funcs = [ 0., 0.8, 0.]
u = np.piecewise(x, conds, funcs)

# u = np.sin(2*np.pi*x/Lx) # sin function only as test
plt.xlim(-Lx, Lx) # x-axis limits acc. to Lx
plt.ylim(np.min(u)-1, np.max(u)+1) # y-axis limits acc. to u

## animation (for n in range(nt))
def animate(n):
    un = u.copy()
    umax = np.max(np.absolute(u))
    dt = dx*Co/umax

    # tn = t
    # t = dt
    # print(tn)

    for i in range(len(x)):
        if umax > 2: break # to avoid explosion
        u[i] = un[i] - un[i]*dt/dx*(un[i] - un[i-1])
    
    line.set_data(x, u)
    frame_text.set_text('frame = %4d \ndt = %.6f' % (n, dt))
    return line, frame_text,

anim = animation.FuncAnimation(fig, animate, frames=1000, interval=10, repeat = False,
                               blit=True) # no init function

## plot setup and run!
plt.title('Step 2: 1-D Non-Linear Convection')
plt.xlabel('x-axis')
plt.ylabel('wave amplitude')
plt.show()

## --------- NOT USED -------------
## testing my function :)
# u = np.exp(-2*x**2) # my function! (as test an exponential growth)
# plt.plot(x, u)

## done!
# plt.show()
# time = datetime.datetime.now()
# print 'done @ %s' % (time.strftime('%H:%M:%S'))

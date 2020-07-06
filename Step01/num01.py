## foreword
# trying again to gain some understanding of the most basic Navie-Stokes cases, August 2019
# based on the animation example from Jake Vanderplas (http://jakevdp.github.com)

## preamble
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

## geometry setup
Lx = 10.0 # length given as float
nx = 1001 # 1D mesh size
dx = Lx/(nx-1) # space step size
Ly = 1 # wave max amplitude + a bit more

## figure setup
fig = plt.figure()
ax = plt.axes(xlim=(0, Lx), ylim=(0, Ly))
line, = ax.plot([], [], lw=2) # why is it a comma after the onject declaration?
frame_text = ax.text(0.75, 0.8, '', transform=ax.transAxes) # for printing text in the plot

## wave setup
c = 1 # wave velocity
Co = 1 # Courant number, should be =<1 to guarantee numerical stability
# t = 0.0 # time
t_tot = 10 # total time in seconds
dt = dx*Co/c # time step acc. to the CLF stability condition, used as "interval" in FuncAnimation
nt = int(t_tot/dt) # total number of iterations, used as "frames" in FuncAnimation
x = np.linspace(0, Lx, nx)
u = np.zeros(nx) # init in 0
# u[int(0.5/dx):int(1.5/dx)] = 2 # "wave" interval with value x
u = np.sin(2*np.pi*x/Lx) # sin function, period = Lx
plt.ylim(np.min(u)-0.1, np.max(u)+0.1) # resize of the y-axis limits 

# ## initialization function (not really used in this example)
# def init():
#     line.set_data([], [])
#     return line, # don't forget the comma!
 
## animation function (called sequentially)
def animate(n): # for n in range(nt)
    # u[int(2.5/dx):int(5.0/dx)] = 0 + (n/100.0) # just testing, but don't forget the floats!
    t =  n*dt
    un = u.copy()
    # for i in range(1,nx): # 0 not included in the loop
    for i in range(len(x)): # 0 included in the loop    
        u[i] = un[i] - c*dt/dx*(un[i]-un[i-1])
        
    line.set_data(x, u)
    frame_text.set_text('frame = %s \ntime = %.6s \ndt = %.6s' % (n, t, dt))
    return line, frame_text # don't forget the comma!

## call animator (if blit=True, only changed parts are redrawn)
# anim = animation.FuncAnimation(fig, animate, init_func=init, frames=nt, interval=20, blit=True)
anim = animation.FuncAnimation(fig, animate, frames=nt, interval=dt*1000, repeat = False,
                               blit=True) # no init function

## plot setup and run!
plt.title('Step 1: 1-D Linear Convection')
# plt.text(1, 0.65, 'frame = %s'%(anim.frame_seq)) # not working!! 
plt.xlabel('x-axis')
plt.ylabel('wave amplitude')
plt.show()

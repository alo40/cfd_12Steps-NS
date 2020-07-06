"""
My user-defined functions that will be called from other files.
More info here: https://problemsolvingwithpython.com/07-Functions-and-Modules/07.05-Calling-Functions-from-Other-Files/
Today's 14.06.2020, 2nd wave of Quarentene is coming ...
"""

# convection time step acc. to the CLF stability criterion
# ------------------------------------------------------
def convection_timeStep(cx, cy, dx, dy, Co):
    
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

    return min(dtx, dty)      # total time step, choose min value to assure numerical stability


# 2D non-linear convection
# ------------------------------------------------------
def convection_stencil(u, v, dx, dy, dt):
    
    # copy u, v
    un = u.copy()
    vn = v.copy()
    
    # 2D array operators (non-code)
    # u(j  , i  ) --> u[1:-1, 1:-1]
    # u(j  , i-1) --> u[1:-1, 0:-2]
    # u(j-1, i  ) --> u[0:-2, 1:-1]
    
    # u stencil
    u[1:-1, 1:-1] = (un[1:-1, 1:-1] - (un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, 0:-2])) -
                                      (vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[0:-2, 1:-1])))

    # v stencil
    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - (un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2])) -
                                      (vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1])))

    # u boundary condtion (cyclic)
    u[ 0,  :] = u[-1,  :]
    u[-1,  :] = u[-2,  :]
    u[ :,  0] = u[ :, -1]
    u[ :, -1] = u[ :, -2]
    
    # v boundary condtion (cyclic)
    v[ 0,  :] = v[-1,  :]
    v[-1,  :] = v[-2,  :]
    v[ :,  0] = v[ :, -1]
    v[ :, -1] = v[ :, -2]
    
    return u, v # simple as that!


# diffusion time step acc. to the CLF stability criterion
# ------------------------------------------------------
def diffusion_timeStep(dx, dy, Co, Nu):

    # x-time step
    dtx = Co * dx**2 / Nu

    # y-time step
    dty = Co * dy**2 / Nu      

    return min(dtx, dty)      # total time step, choose min value to assure numerical stability


# 2D diffusion
# ------------------------------------------------------
def diffusion_stencil(u, dx, dy, dt, Nu):
    
    # copy u
    un = u.copy()

    # 2D array operators (non-code)
    # u(j  , i+1) --> u[1:-1, 2:  ]
    # u(j  , i-1) --> u[1:-1, 0:-2]
    # u(j  , i  ) --> u[1:-1, 1:-1]
    # u(j-1, i  ) --> u[0:-2, 1:-1]
    # u(j+1, i  ) --> u[2:  , 1:-1]

    # u stencil
    u[1:-1, 1:-1] = (un[1:-1, 1:-1] + (Nu * dt / dx**2 * (un[1:-1, 2:  ] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2])) +
                                      (Nu * dt / dy**2 * (un[2:  , 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

    # u boundary condtion (cyclic)
    u[ 0,  :] = u[-1,  :]
    u[-1,  :] = u[-2,  :]
    u[ :,  0] = u[ :, -1]
    u[ :, -1] = u[ :, -2]

    return u

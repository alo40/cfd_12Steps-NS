"""
My user-defined functions that will be called from other files.
More info here: https://problemsolvingwithpython.com/07-Functions-and-Modules/07.05-Calling-Functions-from-Other-Files/
Today's 21.06.2020, need to find a way to have only one file for all my functions ...
"""

# convection time step acc. to the CLF stability criterion
# ------------------------------------------------------
def timeStep_convection(cx, cy, dx, dy, Co):
    
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


# diffusion time step acc. to the CLF stability criterion
# ------------------------------------------------------
def timeStep_diffusion(dx, dy, Co, Nu):

    # x-time step
    dtx = Co * dx**2 / Nu

    # y-time step
    dty = Co * dy**2 / Nu      

    return min(dtx, dty)      # total time step, choose min value to assure numerical stability


# Burgers' equation time step acc. to the CLF stability criterion
# ------------------------------------------------------
def timeStep_Burgers(cx, cy, dx, dy, Co, Nu):
    
    # convection
    dt1 = timeStep_convection(cx, cy, dx, dy, Co)
     
    # diffusion
    dt2 = timeStep_diffusion(dx, dy, Co, Nu)

    return min(dt1, dt2)      # total time step, choose min value to assure numerical stability


# Navier-Stokes equations time step acc. to the CLF stability criterion
# ------------------------------------------------------
def timeStep_NS(cx, cy, dx, dy, Co, Nu, rho):
    
    # convection
    dt1 = timeStep_convection(cx, cy, dx, dy, Co)
     
    # diffusion
    dt2 = timeStep_diffusion(dx, dy, Co, Nu)
    
    # pressure
    dt3 = Co * rho * 2 * dy 
    
    return min(dt1, dt2, dt3) # total time step, choose min value to assure numerical stability


# constant boundary conditions
# ------------------------------------------------------
def boundaries_constant(u, a):

    # u boundary condition 
    u[ 0,  :] = a
    u[-1,  :] = a
    u[ :,  0] = a
    u[ :, -1] = a

    return u # simple as that!


# derivative boundary conditions
# ------------------------------------------------------
def boundaries_derivative(u):

    # u boundary condition 
    u[ 0,  :] = u[ 1,  :]
    u[-1,  :] = u[-2,  :]
    u[ :,  0] = u[ :,  1]
    u[ :, -1] = u[ :, -2]

    return u # simple as that!    


# cyclic boundary conditions
# ------------------------------------------------------
def boundaries_cyclic(u):
    
    # u boundary condition 
    u[ 0,  :] = u[-1,  :]
    u[-1,  :] = u[-2,  :]
    u[ :,  0] = u[ :, -1]
    u[ :, -1] = u[ :, -2]
    
    return u # simple as that!


# Laplace boundary conditions
# ------------------------------------------------------
def boundaries_Sihi(p):

    # @ x = 0  
    p[ 0,  :] = 1
    
    # @ x = Lx
    p[-1,  :] = 1
    
    # @ y = 0
    p[ :,  0] = 1
    
    # @ y = 0
    p[ :, -1] = 1
    
    return p # simple as that!


# Laplace boundary conditions
# ------------------------------------------------------
def boundaries_Laplace(p, y):

    # @ x = 0  
    p[0, :] = y # gradient
    
    # @ x = Lx
    p[-1, :] = 0
    
    # @ y = 0
    p[:, 0] = p[:, 1] # dp/dy = 0
    
    # @ y = 0
    p[:, -1] = p[:, -2] # dp/dy = 0
    
    return p # simple as that!


# 2D non-linear convection
# ------------------------------------------------------
def stencil_convection(u, v, dx, dy, dt):
    
    # copy u, v
    un = u.copy()
    vn = v.copy()
       
    # u stencil
    u[1:-1, 1:-1] = (un[1:-1, 1:-1] - (un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, 0:-2])) -
                                      (vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[0:-2, 1:-1])))

    # v stencil
    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - (un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2])) -
                                      (vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1])))

    return u, v # simple as that!


# 2D diffusion
# ------------------------------------------------------
def stencil_diffusion(u, dx, dy, dt, Nu):
    
    # copy u
    un = u.copy()

    # u stencil
    u[1:-1, 1:-1] = (un[1:-1, 1:-1] + (Nu * dt / dx**2 * (un[1:-1, 2:  ] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2])) +
                                      (Nu * dt / dy**2 * (un[2:  , 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

    return u # simple as that!


# 2D Burgers' equation
# ------------------------------------------------------
def stencil_Burgers(u, v, dx, dy, dt, Nu):
        
    # copy u, v
    un = u.copy()
    vn = v.copy()
    
    # call convection stencil    
    u1, v1 = stencil_convection(u, v, dx, dy, dt)
    
    # call diffusion stencils
    u2 = stencil_diffusion(u, dx, dy, dt, Nu)
    v2 = stencil_diffusion(v, dx, dy, dt, Nu)
    
    # convection + diffusion      
    u = u1 + u2
    v = v1 + v2
    
    # correct extra term u(i,j), v(i,j)
    u[1:-1, 1:-1] = u[1:-1, 1:-1] - un[1:-1, 1:-1] 
    v[1:-1, 1:-1] = v[1:-1, 1:-1] - vn[1:-1, 1:-1]
    
    return u, v # simple as that!


# 2D Laplace equation
# ------------------------------------------------------
def stencil_Laplace(p, dx, dy):
    
    from numpy.linalg import norm
    
    # copy p
    pn = p.copy()
    
    # p stencil
    p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:  ] + pn[1:-1, 0:-2]) + 
                      dx**2 * (pn[2:  , 1:-1] + pn[0:-2, 1:-1])) / (2 * (dx**2 + dy**2)))
    
    # norm difference
    L = norm(p) - norm(pn)
    
    return p, L # simple as that!


# 2D Laplace equation
# ------------------------------------------------------
def stencil_Poisson(p, b, dx, dy):
    
    from numpy.linalg import norm
    
    # copy p
    pn = p.copy()
        
    # p stencil
    p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:  ] + pn[1:-1, 0:-2]) + 
                      dx**2 * (pn[2:  , 1:-1] + pn[0:-2, 1:-1]) -
                      dx**2 * dy**2 * b[1:-1, 1:-1]) / (2 * (dx**2 + dy**2)))
    
    # norm difference
    L = norm(p) - norm(pn)
    
    return p, L # simple as that!


# 2D Navier-Stokes equations without source term
# ------------------------------------------------------
def stencil_NS_noSource(u, v, p, dx, dy, dt, Nu, rho):
    
    # copy all vectors
    un = u.copy()
    vn = v.copy()
    pn = p.copy()
    
    # Calculate pressure using div(u)=0 (Poisson)
    p[1:-1, 1:-1] = (
        #
                    1/(2*(dx**2 + dy**2))
                  * (
        #
                    (pn[1:-1, 2:  ] + pn[1:-1, 0:-2])*dy**2 
                  + (pn[2:  , 1:-1] + pn[0:-2, 1:-1])*dx**2                   
        #
                  - dx**2*dy**2*rho
                  * (
        #
                  - ((un[1:-1, 2:  ] - un[1:-1, 0:-2])/(2*dx))**2
                  - 2*((un[2:  , 1:-1] - un[0:-2, 1:-1])/(2*dy))*((vn[1:-1, 2:  ] - vn[1:-1, 0:-2])/(2*dx))
                  - ((vn[2:  , 1:-1] - vn[0:-2, 1:-1])/(2*dy))**2
        #
                  + 1/dt
                  * (
        #
                    (un[1:-1, 2:  ] - un[1:-1, 0:-2])/(2*dx)
                  + (vn[2:  , 1:-1] - vn[0:-2, 1:-1])/(2*dy) 
        #
                    ))))
    
    # Correct velocity using pressure
    u[1:-1, 1:-1] = (
        #
                    un[1:-1, 1:-1]
                  - un[1:-1, 1:-1]*dt/dx*(un[1:-1, 1:-1] - un[1:-1, 0:-2])
                  - vn[1:-1, 1:-1]*dt/dy*(un[1:-1, 1:-1] - un[0:-2, 1:-1])
        #
                  - dt/(2*dx*rho)*(p[1:-1, 2:  ] - p[1:-1, 0:-2]) # pressure using div(u)=0  
        #
                  + Nu*dt/dx**2*(un[1:-1, 2:  ] - 2*un[1:-1, 1:-1] + un[1:-1, 0:-2])
                  + Nu*dt/dy**2*(un[2:  , 1:-1] - 2*un[1:-1, 1:-1] + un[0:-2, 1:-1])
        #
                    )
    
    v[1:-1, 1:-1] = (
        #
                    vn[1:-1, 1:-1]
                  - un[1:-1, 1:-1]*dt/dx*(vn[1:-1, 1:-1] - vn[1:-1, 0:-2])
                  - vn[1:-1, 1:-1]*dt/dy*(vn[1:-1, 1:-1] - vn[0:-2, 1:-1])
        #
                  - dt/(2*dy*rho)*(p[2:  , 1:-1] - p[0:-2, 1:-1])
        #
                 + Nu*dt/dx**2*(vn[1:-1, 2:  ] - 2*vn[1:-1, 1:-1] + vn[1:-1, 0:-2])
                 + Nu*dt/dy**2*(vn[2:  , 1:-1] - 2*vn[1:-1, 1:-1] + vn[0:-2, 1:-1])
        #
                    )
    
    return u, v, p # not so simple this time!
    

# 2D array operators (non-code)
# ------------------------------------------------------
# u(j  , i+1) --> u[1:-1, 2:  ]
# u(j  , i-1) --> u[1:-1, 0:-2]
# u(j  , i  ) --> u[1:-1, 1:-1]
# u(j-1, i  ) --> u[0:-2, 1:-1]
# u(j+1, i  ) --> u[2:  , 1:-1]                    
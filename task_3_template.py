# MTF072 Computational Fluid Dynamics
# Task 3: laminar lid-driven cavity
# Template prepared by:
# Gonzalo Montero Villar
# Department of Mechanics and Maritime Sciences
# Division of Fluid Dynamics
# villar@chalmers.se
# December 2019

#==============Packages needed=================
import matplotlib.pyplot as plt
import numpy as np

#================= Inputs =====================

# Fluid properties and B. C. inputs

UWall =  # velocity of the upper wall
rho   =  # density
nu    =  # kinematic viscosity

data_file = # data file where the given solution is stored

# Geometric inputs (fixed so that a fair comparison can be made)

nI = 11 # number of nodes X direction. 
nJ = 11 # number of nodes Y direction. 
xL =  1 # length in X direction
yL =  1 # length in Y direction

# Solver inputs

nIterations           =  # maximum number of iterations
n_inner_iterations_gs =  # amount of inner iterations when solving 
                              # pressure correction with Gauss-Seidel
resTolerance =  # convergence criteria for residuals
                     # each variable
alphaUV =       # under relaxation factor for U and V
alphaP  =       # under relaxation factor for P

# ================ Code =======================

# For all the matrices the first input makes reference to the x coordinate
# and the second input to the y coordinate, (i+1) is east and (j+1) north

# Allocate all needed variables
cI = nI + 1                      # number of cells in the X direction. Cells 
                                  # added in the boundaries
cJ = nJ + 1                      # number of cells in the Y direction. Cells 
                                  # added in the boundaries
coeffsUV   = np.zeros((cI,cJ,5)) # coefficients for the U and V equation
                                  # E, W, N, S and P
sourceUV   = np.zeros((cI,cJ,2)) # source coefficients for the U and V equation
                                  # U and V
coeffsPp   = np.zeros((cI,cJ,5)) # coefficients for the pressure correction
                                  # equation E, W, N, S and P
sourcePp   = np.zeros((cI,cJ))   # source coefficients for the pressure
                                  # correction equation
U          = np.zeros((cI,cJ))   # U velocity matrix
V          = np.zeros((cI,cJ))   # V velocity matrix
P          = np.zeros((cI,cJ))   # pressure matrix
Pp         = np.zeros((cI,cJ))   # pressure correction matrix

massFlows  = np.zeros((cI,cJ,4)) # mass flows at the faces
                                  # m_e, m_w, m_n and m_s

residuals  = np.zeros((3,1))     # U, V and conitnuity residuals

# Generate mesh and compute geometric variables

# Allocate all variables matrices
xCoords_C = np.zeros((cI,cJ)) # X coords of the cells
yCoords_C = np.zeros((cI,cJ)) # Y coords of the cells
xCoords_N = np.zeros((nI,nJ)) # X coords of the nodes
yCoords_N = np.zeros((nI,nJ)) # Y coords of the nodes
dxe_C     = np.zeros((cI,cJ)) # X distance to east cell
dxw_C     = np.zeros((cI,cJ)) # X distance to west cell
dyn_C     = np.zeros((cI,cJ)) # Y distance to north cell
dys_C     = np.zeros((cI,cJ)) # Y distance to south cell
dx_C      = np.zeros((cI,cJ)) # X size of the cell
dy_C      = np.zeros((cI,cJ)) # Y size of the cell

residuals_U = []
residuals_V = []
residuals_c = []

dx = xL/(nI - 1)
dy = yL/(nJ - 1)

# Fill the coordinates
for i in range(nI):
    for j in range(nJ):
        # For the nodes
        xCoords_N[i,j] = i*dx
        yCoords_N[i,j] = j*dy

        # For the cells
        if i > 0:
            xCoords_C[i,j] = 0.5*(xCoords_N[i,j] + xCoords_N[i-1,j])
        if i == nI-1 and j>0:
            yCoords_C[i+1,j] = 0.5*(yCoords_N[i,j] + yCoords_N[i,j-1])
        if j > 0:
            yCoords_C[i,j] = 0.5*(yCoords_N[i,j] + yCoords_N[i,j-1])
        if j == nJ-1 and i>0:
            xCoords_C[i,j+1] = 0.5*(xCoords_N[i,j] + xCoords_N[i-1,j])

        # Fill dx_C and dy_C
        if i > 0:
            dx_C[i,j] = xCoords_N[i,j] - xCoords_N[i-1,j]
        if j > 0:
            dy_C[i,j] = yCoords_N[i,j] - yCoords_N[i,j-1]

xCoords_C[-1,:] = xL
yCoords_C[:,-1] = yL


# Fill dxe, dxw, dyn and dys
for i in range(1,cI-1):
    for j in range(1,cJ-1):
        dxe_C[i,j] = xCoords_C[i+1,j] - xCoords_C[i,j]
        dxw_C[i,j] = xCoords_C[i,j] - xCoords_C[i-1,j]
        dyn_C[i,j] = yCoords_C[i,j+1] - yCoords_C[i,j]
        dys_C[i,j] = yCoords_C[i,j] - yCoords_C[i,j-1]

# Initialize variable matrices

U[:,:] = 
V[:,:] = 
P[:,:] = 

# Looping

for iter in range(nIterations):
    # Impose boundary conditions for velocities, only the top boundary wall
    # is moving from left to right with UWall
    
    # Impose pressure boundary condition, all homogeneous Neumann
    
    # Compute coefficients for U and V equations
    for i in range(1,cI-1):
        for j in range(1,cJ-1):
            coeffsUV[i,j,0] = 
            coeffsUV[i,j,1] = 
            coeffsUV[i,j,2] = 
            coeffsUV[i,j,3] = 
            coeffsUV[i,j,4] = 
            sourceUV[i,j,0] = 
            sourceUV[i,j,1] = 
    
    # Introduce implicit under-relaxation for U and V
    for i in range(1,cI-1):
        for j in range(1,cJ-1):
  
    # Solve for U and V using Gauss-Seidel

    # Calculate at the faces using Rhie-Chow for the face velocities
    
    # Calculate pressure correction equation coefficients
    
    for i in range(1,cI-1):
        for j in range(1,cJ-1):
            
            coeffsPp[i,j,0] = 
            coeffsPp[i,j,1] = 
            coeffsPp[i,j,2] = 
            coeffsPp[i,j,3] = 
            coeffsPp[i,j,4] = 
            sourcePp[i,j]   = 
    
    # Solve for pressure correction (Note that more that one loop is used)
    for iter_gs in range(n_inner_iterations_gs):
        for j in range(1,cJ-1):
            for i in range(1,cI-1):    
    
    # Set Pp with reference to cell (2,2) and copy to boundaries
    
    # Correct velocities, pressure and mass flows

    for i in range(1,cI-1):
        for j in range(1,cJ-1):         
    
    # impose zero mass flow at the boundaries

    # Copy P to boundaries
    
    # Compute residuals
    residuals_U.append(0) # U momentum residual
    residuals_V.append(0) # V momentum residual
    residuals_c.append(0) # continuity residual

    for i in range(1,cI-1):
        for j in range(1,cJ-1):
            residuals_U[-1] = 
            residuals_V[-1] = 
            residuals_c[-1] = 

    print('iteration: %d\nresU = %.5e, resV = %.5e, resCon = %.5e\n\n'\
        % (iter, residuals_U[-1], residuals_V[-1], residuals_c[-1]))
    
    #  Check convergence
    if resTolerance>max([residuals_U[-1], residuals_V[-1], residuals_c[-1]]):
        break

# Plotting section (these are some examples, more plots might be needed)


# Plot mesh
plt.figure()
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Computational mesh')

# Plot results

plt.figure()

# U velocity contour
plt.subplot(2,3,1)
plt.title('U velocity [m/s]')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

# V velocity contour
plt.subplot(2,3,2)
plt.title('V velocity [m/s]')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

# P contour
plt.subplot(2,3,3)
plt.title('Pressure [Pa]')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

# Vector plot
plt.subplot(2,3,4)
plt.title('Vector plot of the velocity field')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

# Comparison with data
data=np.genfromtxt(data_file, skip_header=1)
uInterp = np.zeros((cJ-2,1))
vInterp = np.zeros((cJ-2,1))
for j in range(1,cJ-1):
    for i in range(1,cI-1):
        if xCoords_C[i,j]<0.5 and xCoords_C[i+1,j]>0.5:
            uInterp[j-1] = (U[i+1,j] + U[i,j])*0.5
            vInterp[j-1] = (V[i+1,j] + V[i,j])*0.5
            break
        elif abs(xCoords_C[i,j]-0.5) < 0.000001:
            uInterp[j-1] = U[i,j]
            vInterp[j-1] = V[i,j]
            break

plt.subplot(2,3,5)
plt.plot(data[:,0],data[:,2],'r.',markersize=20,label='data U')
plt.plot(data[:,1],data[:,2],'b.',markersize=20,label='data V')
plt.plot(uInterp,yCoords_C[1,1:-1],'k',label='sol U')
plt.plot(vInterp,yCoords_C[1,1:-1],'g',label='sol V')
plt.title('Comparison with data at x = 0.5')
plt.xlabel('u, v [m/s]')
plt.ylabel('y [m]')
plt.legend()

plt.subplot(2,3,6)
plt.title('Residual convergence')
plt.xlabel('iterations')
plt.ylabel('residuals [-]')
#plt.legend('U momentum','V momentum', 'Continuity')
plt.title('Residuals')
plt.show()




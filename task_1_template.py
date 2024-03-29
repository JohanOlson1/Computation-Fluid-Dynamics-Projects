# MTF073 Computational Fluid Dynamics
# Task 1: 2D diffusion equation
# Template prepared by:
# Gonzalo Montero Villar
# Department of Mechanics and Maritime Sciences
# Division of Fluid Dynamics
# villar@chalmers.se
# November 2020
# Packages needed
import numpy as np
import matplotlib.pyplot as plt

#===================== Schematic ==================
#
#                  0----------------0
#                  |                |
#                  |                |
#                  |    [i,j+1]     |
#                  |       X        |
#                  |      (N)       |
#                  |                |
#                  |                |
# 0----------------0----------------0----------------0
# |                |                |                |
# |                |                |                |
# |    [i-1,j]     |     [i,j]      |    [i+1,j]     |
# |       X        |       X        |       X        |
# |      (W)       |      (P)       |      (E)       |
# |                |                |                |
# |                |                |                |
# 0----------------0----------------0----------------0
#                  |                |
#                  |                |
#                  |    [i,j-1]     |
#                  |       X        |
#                  |      (S)       |
#                  |                |
#                  |                |
#                  0----------------0
#
# X:  marks the position of the nodes, which are the centers
#     of the control volumes, where temperature is computed.
# 0:  marks the position of the mesh points or control volume
#     corners.
# []: in between square brakets the indexes used to find a 
#     node with respect to the node "P" are displayed.
# (): in between brakets the name given to refer to the nodes
#     in the lectures as well as in the book with respect to 
#     the node "P" are displayed.   

#===================== Inputs =====================

# Constants

T1 = 10
T2 = 20
c1 = 25
c2 = 0.25
 
# Geometric inputs

mI = 4 # number of mesh points X direction.
mJ = 4 # number of mesh points Y direction.
grid_type = 'equidistant' # this sets equidistant mesh sizing or non-equidistant
xL = 1 # length of the domain in X direction
yL = 1 # length of the domain in Y direction

# Solver inputs

nIterations  = 1 # maximum number of iterations
resTolerance =  0.001 # convergence criteria for residuals each variable

#====================== Code ======================

# For all the matrices the first input makes reference to the x coordinate
# and the second input to the y coordinate (check Schematix above)

# Allocate all needed variables
nI = mI + 1                    # number of nodes in the X direction. Nodes 
                               # added in the boundaries
nJ = mJ + 1                    # number of nodes in the Y direction. Nodes 
                               # added in the boundaries
coeffsT = np.zeros((nI,nJ,5))  # coefficients for temperature
                               # E, W, N, S and P, the a:s
                               
S_U     = np.zeros((nI,nJ))    # source term for temperature
S_P     = np.zeros((nI,nJ))    # source term for temperature
T       = np.zeros((nI,nJ))    # temperature matrix
k       = np.zeros((nI,nJ))    # coefficient of conductivity
q       = np.zeros((nI,nJ,2))  # heat flux, first x and then y component

residuals = [] # List containing the value of the residual for each iteration

# Generate mesh and compute geometric variables

# Allocate all variables matrices
xCoords_M = np.zeros((mI,mJ)) # X coords of the mesh points
yCoords_M = np.zeros((mI,mJ)) # Y coords of the mesh points
xCoords_N = np.zeros((nI,nJ)) # X coords of the nodes
yCoords_N = np.zeros((nI,nJ)) # Y coords of the nodes
dxe_N     = np.zeros((nI,nJ)) # X distance to east node
dxw_N     = np.zeros((nI,nJ)) # X distance to west node
dyn_N     = np.zeros((nI,nJ)) # Y distance to north node
dys_N     = np.zeros((nI,nJ)) # Y distance to south node
dx_CV     = np.zeros((nI,nJ)) # X size of the control volume
dy_CV     = np.zeros((nI,nJ)) # Y size of the control volume

if grid_type == 'equidistant':
    # Control volume size
    dx = xL/(mI - 1)
    dy = yL/(mJ - 1)

    # Fill the coordinates
    for i in range(mI):
        for j in range(mJ):
            # For the mesh points
            xCoords_M[i,j] = i*dx
            yCoords_M[i,j] = j*dy

            # For the nodes
            if i > 0:
                xCoords_N[i,j] = 0.5*(xCoords_M[i,j] + xCoords_M[i-1,j]) # All x node coordinates
            if i == (mI-1) and j>0:
                yCoords_N[i+1,j] = 0.5*(yCoords_M[i,j] + yCoords_M[i,j-1]) # Last row of y
            if j >0:
                yCoords_N[i,j] = 0.5*(yCoords_M[i,j] + yCoords_M[i,j-1]) # All y node coord.
            if j == (mJ-1) and i>0:
                xCoords_N[i,j+1] = 0.5*(xCoords_M[i,j] + xCoords_M[i-1,j]) # last column of x, we miss without this one

            # Fill dx_CV and dy_CV
            if i>0:
                dx_CV[i,j] = xCoords_M[i,j] - xCoords_M[i-1,j]
            if j>0:
                dy_CV[i,j] = yCoords_M[i,j] - yCoords_M[i,j-1]
                

# =============================================================================
# elif grid_type == 'non-equidistant':
#     rx = 1.15
#     ry = 1.15
#     
#     # Fill the necessary code to generate a non equidistant grid and
#     # fill the needed matrixes for the geometrical quantities
#     
# xCoords_N[-1,:] = xL
# yCoords_N[:,-1] = yL
# 
# #TODO: FIX
# # Fill dxe, dxw, dyn and dys
# for i in range(1,nI - 1):
#     for j in range(1,nJ - 1):
#         dxe_N[i,j] = 
#         dxw_N[i,j] = 
#         dyn_N[i,j] = 
#         dys_N[i,j] = 
# =============================================================================


# Initialize variable matrices and boundary conditions
                
    for i in range(2,nI-2):
        T[i,0] = T1                                   #B1
        T[i,nJ-1] = 5 + 3*(1+5*xCoords_N[i,nJ-1]/xL)  #B3

    for j in range(2,nJ-2):
        T[nI-1,j] = T2                                #B2
        # What about B4?!

# Looping

for iter in range(nIterations):
    
    # Update conductivity coefficient matrix, k, according to your case
    for i in range(nI):
        for j in range(NJ):
            k[i,j] = 2(1+20*T[i,j]/T1)
        
            # Update source term matrix according to your case
            S_U = 15(c1-c2*T[i,j]**2) * dx * dy    # In loop
    
    # Compute coeffsT for all the nodes which are not boundary nodes
    ## Compute coefficients for nodes one step inside the domain
    
    ### First, north and south boundaries
    for i in range(2,nI-2):
        # Why start at 2?! Should be 1?
        coeffsT[i,1,0] =  

    ### Second, east and west boundaries,  2, 4
    for j in range(2,nJ-2):
 
        
    ## Compute coefficients for inner nodes
    #TODO: FIX
    for i in range(2,nI-2):
        for j in range(2,nJ-2): # This is not done for all!!!!!!, we also need the ones bordering the boundary
            coeffsT[i,j,0] = k_w * (yCoords_M[i,j] - yCoords_M[i,j-1])/(xCoords_N[i,j-1] - xCoords_N[i-1,j-1])
            coeffsT[i,j,1] = k_e * (yCoords_M[i+1,j] - yCoords_M[i+1,j-1])/(xCoords_N[i+1,j-1] - xCoords_N[i,j-1])
            coeffsT[i,j,2] = k_s * (xCoords_M[i,j] - xCoords_M[i-1,j])/(yCoords_N[i,j] - yCoords_N[i,j-1])
            coeffsT[i,j,3] = k_n * (xCoords_M[i,j+1] - xCoords_M[i-1,j+1])/(yCoords_N[i,j+1] - yCoords_N[i,j])
            coeffsT[i,j,4] = coeffsT[i,j,0] + coeffsT[i,j,1] + coeffsT[i,j,2] + coeffsT[i,j,3] 
    
    ## Compute coefficients corner nodes (one step inside)

    # Solve for T using Gauss-Seidel
    
    # Copy T to boundaries where homegeneous Neumann needs to be applied
    
    # Compute residuals (taking into account normalization)
    r = 0
    
    residuals.append(r)
    
    print('iteration: %d\nresT = %.5e\n\n'  % (iter, residuals[-1]))
    
    #  Check convergence
    if resTolerance>residuals[-1]:
        break

# Compute heat fluxes
for i in range(1,nI-1):
    for j in range(1,nJ-1):
        q[i,j,0] = 
        q[i,j,1] = 
    
# Plotting section (these are some examples, more plots might be needed)

# Plot results
plt.figure()

# Plot mesh
plt.subplot(2,2,1)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Computational mesh')
plt.axis('equal')

# Plot temperature contour
plt.subplot(2,2,2)
plt.title('Temperature [ºC]')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')

# Plot residual convergence
plt.subplot(2,2,3)
plt.title('Residual convergence')
plt.xlabel('iterations')
plt.ylabel('residuals [-]')
plt.title('Residual')

# Plot heat fluxes
plt.subplot(2,2,4)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Heat flux')
plt.axis('equal')

plt.show()

    



# MTF072 Computational Fluid Dynamics
# Task 2: convection-diffusion equation
# Template prepared by:
# Gonzalo Montero Villar
# Department  of Mechanics and Maritime Sciences
# Division of Fluid Dynamics
# villar@chalmers.se
# November 2019

# The script assumes that the folder with data is in the same path as this file

# Packages needed
import numpy as np
import matplotlib.pyplot as plt

def Discretize(Pe,):
    if Pe < 2:
        # central
    else:
        # upwind
    return

def  ReadDataAndGeometry(caseID, grid_type):
	if caseID <= 5:
	    grid_number = 1
	elif caseID <= 10:
	    grid_number = 2
	elif caseID <= 15:
	    grid_number = 3
	elif caseID <= 20:
	    grid_number = 4
	elif caseID <= 25:
	    grid_number = 5

	path = 'data/grid%d/%s_grid' % (grid_number,grid_type)

	# Read data
	xCoords_M = np.genfromtxt('%s/xc.dat' % (path)) # x node coordinates
	yCoords_M = np.genfromtxt('%s/yc.dat' % (path)) # y node coordinates
	u_data = np.genfromtxt('%s/u.dat' % (path))     # U velocity at the nodes
	v_data = np.genfromtxt('%s/v.dat' % (path))     # V veloctiy at the nodes

	# Allocate geometrical data and variables
	mI        = len(xCoords_M);      # number of mesh points in the x direction
	mJ        = len(yCoords_M);      # number of mesh points in the y direction
	cI        = mI + 1;              # number of nodes in the x direction
	cJ        = mJ + 1;              # number of nodes in the y direction
	xCoords_N = np.zeros((cI,1));    # mesh points x coordinates
	yCoords_N = np.zeros((cJ,1));    # mesh points y coordinates
	dxe_C     = np.zeros((cI,1));    # X distance to east cell
	dxw_C     = np.zeros((cI,1));    # X distance to west cell
	dyn_C     = np.zeros((cJ,1));    # Y distance to north cell
	dys_C     = np.zeros((cJ,1));    # Y distance to south cell
	dx_C      = np.zeros((cI,1));    # X size of the cell
	dy_C      = np.zeros((cJ,1));    # Y size of the cell

	# Fill the cell coordinates as the mid point between mesh points, and at the same
	# position at the boundaries. Compute cell sizes
	for i in range(1,cI-1):
	    xCoords_N[i] = (xCoords_M[i] + xCoords_M[i-1])/2
	    dx_C[i]      = xCoords_M[i] - xCoords_M[i-1]
	
	for j in range(1,cJ-1):
	    yCoords_N[j] = (yCoords_M[j] + yCoords_M[j-1])/2
	    dy_C[j]      = yCoords_M[j] - yCoords_M[j-1]
	
	xCoords_N[0]  = xCoords_M[0]
	xCoords_N[-1] = xCoords_M[-1]
	yCoords_N[0]  = yCoords_M[0]
	yCoords_N[-1] = yCoords_M[-1]

	# Compute distances between nodes
	for i in range(1,cI-1): 
	    dxe_C[i] = xCoords_N[i+1] - xCoords_N[i]
	    dxw_C[i] = xCoords_N[i] - xCoords_N[i-1]
	
	for j in range(1,cJ-1):
	    dyn_C[j] = yCoords_N[j+1] - yCoords_N[j]
	    dys_C[j] = yCoords_N[j] - yCoords_N[j-1]
	

	# Reshape the velocity data
	U = u_data.reshape(cI,cJ)
	V = v_data.reshape(cI,cJ)

	return [xCoords_M, yCoords_M, mI, mJ, cI, cJ, xCoords_N, yCoords_N, dxe_C, dxw_C, dyn_C, dys_C, dx_C, dy_C, U, V]

# Inputs

grid_type   = 'coarse' # either 'coarse' or 'fine'
caseID      =  15      # your case number to solve
k           =  1      
rho         =  1      # density
nIterations =  2000     # number of iterations
Cp          =  200
plotVelocityVectors = True
resTolerance = 0.001

# Read data for velocity fields and geometrical quantities

# For all the matrices the first input makes reference to the x coordinate
# and the second input to the y coordinate, (i+1) is east and (j+1) north

[xCoords_M, yCoords_M, mI, mJ, cI, cJ, xCoords_N, yCoords_N, dxe_C, dxw_C, dyn_C, dys_C, dx_C, dy_C, U, V] = ReadDataAndGeometry(caseID, grid_type)

# Plot velocity vectors if required
if plotVelocityVectors:
	plt.figure()
	plt.quiver(xCoords_N, yCoords_N, U.T, V.T)
	plt.title('Velocity vectors')
	plt.xlabel('x [m]')
	plt.ylabel('y [m]')
	plt.show()

# Allocate needed vairables
T = np.zeros((cI, cJ))        # temperature matrix
D = np.zeros((cI, cJ,4))      # diffusive coefficients e, w, n and s
F = np.zeros((cI, cJ,4))      # convective coefficients e, w, n and s
coeffsT = np.zeros((cI,cJ,5)) # hybrid scheme coefficients E, W, N, S, P
Pe = np.zeros(4)

residuals = []

# Code

gamma = k/Cp

## Diffusive and convective coefficient calculations
for i in range(1,cI-1):
	for j in range(1,cJ-1):
        #TODO: fix function to interpolate gamma
	    D[i,j,0] = gamma_e /dx_C[i,j] # east diffusive
	    D[i,j,1] = gamma_w /dxw_C[i,j] # west diffusive
	    D[i,j,2] = gamma_n /dyn_C[i,j] # north diffusive
	    D[i,j,3] = gamma_s /dys_C[i,j] # south diffusive
	    
        # Remark, unsure of indices, U,V in node or mesh?
	    F[i,j,0] = rho * U[i,j] # east convective
	    F[i,j,1] = rho * U[i-1,j] # west convective
	    F[i,j,2] = rho * V[i,j] # north convective
	    F[i,j,3] = rho * V[i,j-1] # south convective

# Hybrid scheme coefficients calculations (taking into account boundary conditions)
for i in range(1,cI-1):
	for j in range(1,cJ-1):
        Pe[0] = F[i,j,0]/D[i,j,0]
        Pe[1] = F[i,j,1]/D[i,j,1]
        Pe[2] = F[i,j,2]/D[i,j,2]
        Pe[3] = F[i,j,3]/D[i,j,3]
        # TODO: if statement, either central or upwind
		coeffsT[i,j,0] = 
		coeffsT[i,j,1] = 
		coeffsT[i,j,2] =
		coeffsT[i,j,3] =
		coeffsT[i,j,4] = coeffsT[i,j,0] + coeffsT[i,j,1] + coeffsT[i,j,2] + coeffsT[i,j,3]

for iter in range(nIterations): 
    # Impose boundary conditions
    
    # Solve for T using Gauss-Seidel or TDMA (both results need to be 
    # presented)
    
    # TDMA
    # TODO: Calculate all A, C' in for loop
    # TODO: Backward differencing in two for loop, one east-wast, one north-south
    
    # Gauss-Seidel 
    # TODO: one for loop, T_p = .... 

    # Copy temperatures to boundaries
    
    # Compute residuals (taking into account normalization)
    residuals.append() # fill with your residual value for the 
                       # current iteration
    
    print('iteration: %d\nresT = %.5e\n\n' % (iter, residuals[-1]))
    
    # Check convergence
    
    if resTolerance>residuals[-1]:
        break


# Plotting (these are some examples, more plots might be needed)
xv, yv = np.meshgrid(xCoords_N, yCoords_N)

plt.figure()
plt.quiver(xv, yv, U.T, V.T)
plt.title('Velocity vectors')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()

plt.figure()
plt.contourf(xv, yv, T.T)
plt.colorbar()
plt.title('Temperature')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()

## TASK 2 TODO
# 1. A convection term added => new equation 
# 2. New type of discretization: Hybrid scheme, Pe > 2 -> Upwind 
# 3. Use of either Gauss-Seidel or TDMA
# 4. Temperature in Kelvin!
# 5. Implement different termination criteria
#
#
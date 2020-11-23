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

# Group 15 Johan Olson, Alexander Rodin
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

	return [xCoords_M, yCoords_M, mI, mJ, cI, cJ, xCoords_N, yCoords_N, \
         dxe_C, dxw_C, dyn_C, dys_C, dx_C, dy_C, U, V]


# Inputs
grid_type   = 'fine'   # either 'coarse' or 'fine'
solver_type = 'TDMA'   # Either 'TDMA' or 'Gauss-Seidel'
caseID      =  15      # your case number to solve
k           =  1       # Conductivity
rho         =  1       # density
nIterations =  1000    # number of iterations
Cp          =  200     # Heat capacity
plotVelocityVectors = False
resTolerance = 0.001   # Error Tolerarance 

# Read data for velocity fields and geometrical quantities
[xCoords_M, yCoords_M, mI, mJ, cI, cJ, xCoords_N, yCoords_N, dxe_C, \
 dxw_C, dyn_C, dys_C, dx_C, dy_C, U, V] = ReadDataAndGeometry(caseID, grid_type)

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

P_x = np.zeros((cI))          # Matrix for TDMA x-wise
P_y = np.zeros((cJ))          # Matrix for TDMA y-wise

Q_x = np.zeros((cI))          # Matrix for TDMA x-wise
Q_y = np.zeros((cJ))          # Matrix for TDMA y-wise

d_x   = np.zeros((cI))        # Matrix for TDMAs Heat generation/neighbouring 
d_y   = np.zeros((cJ))        # Matrix for TDMAs Heat generation/neighbouring

S_U = np.zeros((cI,cJ))       # source term for temperature
T_D = 273.15 + 10             # Temperature at Inlet

residuals = []
nIter = []

# Code
gamma = k/Cp                  # 
q_s = 100                     # Heat flow South Boundary   
Speed_factor = 1
T_New = 273.15 + 100

## Diffusive and convective coefficient calculations
# Remark: No need to keep inside iteration lopp since independent of temp!
for i in range(1,cI-1):
	for j in range(1,cJ-1):
	    D[i,j,0] = gamma * dy_C[j] / dxe_C[i] # east diffusive
	    D[i,j,1] = gamma * dy_C[j] / dxw_C[i] # west diffusive
	    D[i,j,2] = gamma * dx_C[i] / dyn_C[j] # north diffusive
	    D[i,j,3] = gamma * dx_C[i] / dys_C[j] # south diffusive
	    
        # Remark, unsure of indices, U,V in node or mesh?
	    F[i,j,0] = rho * dy_C[j] * U[i,j]/Speed_factor    # east convective
	    F[i,j,1] = rho * dy_C[j] * U[i,j]/Speed_factor    # west convective
	    F[i,j,2] = rho * dx_C[i] * V[i,j]/Speed_factor    # north convective
	    F[i,j,3] = rho * dx_C[i] * V[i,j]/Speed_factor    # south convective
        
# Hybrid scheme coefficients calculations (taking into account boundary conditions)
for i in range(1,cI-1):
	for j in range(1,cJ-1):
		coeffsT[i,j,0] = max([-F[i,j,0], (D[i,j,0]-F[i,j,0])/2, 0])
		coeffsT[i,j,1] = max([F[i,j,1], (D[i,j,1]+F[i,j,1])/2, 0])
		coeffsT[i,j,2] = max([-F[i,j,2], (D[i,j,2]-F[i,j,2])/2, 0])
		coeffsT[i,j,3] = max([F[i,j,3], (D[i,j,3]+F[i,j,3])/2, 0])
		coeffsT[i,j,4] = coeffsT[i,j,0] + coeffsT[i,j,1] + coeffsT[i,j,2] + \
        coeffsT[i,j,3] + F[i,j,0] - F[i,j,1] + F[i,j,2] - F[i,j,3] 
        # S_P is the coefficient that should otherwise be zero, and therefore already done!

## East and West boundary
#Remark: Neumann Cond. -> coefficient is zero
for j in range(1, cJ-1):
    coeffsT[-2,j,0] = 0
    coeffsT[1,j,1] = 0
        
    # Recalculate a_p
    coeffsT[-2,j,4] = coeffsT[-2,j,0] + coeffsT[-2,j,1] + coeffsT[-2,j,2] + \
    coeffsT[-2,j,3] + F[-2,j,0] - F[-2,j,1] + F[-2,j,2] - F[-2,j,3]
        
    # Recalculate a_p
    #S_U[1,j] = max([-F[1,j,1], (D[1,j,1]-F[1,j,1])/2, 0]) * T_New # If dirichlet
    coeffsT[1,j,4] = coeffsT[1,j,0] + coeffsT[1,j,1] + coeffsT[1,j,2] + \
    coeffsT[1,j,3] + F[1,j,0] - F[1,j,1] + F[1,j,2] - F[1,j,3]

# North and (South Boundary constant heat flow)
for i in range(1, cI-1):
    coeffsT[i,1,3] = 0
    coeffsT[i,-2,2] = 0
    S_U[i,1] = q_s * dx_C[i]
    S_U[i,-2] = max([-F[i,-2,2], (D[i,-2,2]-F[i,-2,2])/2, 0]) * T_D
    
for iter in range(nIterations): 
    # Impose boundary conditions? No need?
    
    # TDMA Solution
    if solver_type == 'TDMA':
        # Calculate all P, Q in for loop
        for j in range(1, cJ-1):
            for i in range(1, cI-1):
                d_x[i] = coeffsT[i,j,2] * T[i,j+1] + coeffsT[i,j,3] * T[i,j-1] + S_U[i,j]
                P_x[i] = coeffsT[i,j,0]/(coeffsT[i,j,4] - coeffsT[i,j,1]*P_x[i-1])
                Q_x[i] = (d_x[i] + coeffsT[i,j,1] * Q_x[i-1])/(coeffsT[i,j,4] - coeffsT[i,j,1]*P_x[i-1])
                
            # Backward differencing: East-West
            for k in range(cI-2, 0, -1):
                T[k,j] = P_x[k] * T[k+1,j] + Q_x[k]
                
        for i in range(1, cI-1):
            for j in range(1, cJ-1):
                d_y[j] = coeffsT[i,j,0] * T[i+1,j] + coeffsT[i,j,1] * T[i-1,j] + S_U[i,j]
                P_y[j] = coeffsT[i,j,2]/(coeffsT[i,j,4] - coeffsT[i,j,3]*P_y[j-1])
                Q_y[j] = (d_y[j] + coeffsT[i,j,3] * Q_y[j-1])/(coeffsT[i,j,4] - coeffsT[i,j,3]*P_y[j-1])
                
            # Backward differencing: North-South
            for k in range(cJ-2, 0, -1):
                T[i,k] = P_y[k] * T[i,k+1] + Q_y[k]
            
    # Gauss-Seidel Solution
    if solver_type == 'Gauss-Seidel':
        for j in range(1,cJ-1):
            for i in range(1,cI-1):
                rhs = S_U[i,j] + coeffsT[i,j,0]*T[i+1,j] + coeffsT[i,j,1]*T[i-1,j] + \
                coeffsT[i,j,2]*T[i,j+1] + coeffsT[i,j,3]*T[i,j-1]
                T[i,j] = (rhs)/(coeffsT[i,j,4])

    # Copy temperatures to boundaries
    for i in range(cI):
        T[i,-1] = T_D
        T[i,0] = T[i,1]
    
    for j in range(1, cJ-1):
        T[0,j]  = T[1,j] # If Neumann here
        #T[0,j]  = T_New # If Dirichlet here
        T[-1,j] = T[-2,j]
    
    # Compute residuals (taking into account normalization)
    r = 0
    for j in range(1,cJ-1):
        for i in range(1,cI-1):
            rhs = S_U[i,j] + coeffsT[i,j,0]*T[i+1,j] + coeffsT[i,j,1]*T[i-1,j] + \
            coeffsT[i,j,2]*T[i,j+1] + coeffsT[i,j,3]*T[i,j-1]
            r += np.abs(coeffsT[i,j,4]*T[i,j] - rhs)
            
    
    Inlet = 0
    Outlet = 0
    for i in range(1,cI-1):
        # North Boundary
        Inlet += rho * V[i,-2] * dy_C[-2] * T[i,-2]
    for j in range(1,cJ-1):
        Outlet += rho * U[1,j] * dx_C[1] * T[1,j] + rho * U[-2,j] * dx_C[-2] * T[-2,j]
    
    r /= np.abs(Inlet - Outlet) # r/F
            
    residuals.append(r) # fill with your residual value for the 
    nIter.append(iter)
                       # current iteration
    
    print('iteration: %d\nresT = %.5e\n\n' % (iter, residuals[-1]))
    
    # Check convergence
    
    if resTolerance>residuals[-1]:
        break

# Total Conservation
Flux = 0
for i in range(1,cI-1):
    Flux += - k * (T[i,0] - T[i,1])/dx_C[i]
    Flux += - k * (T[i,-1] - T[i,-2])/dx_C[i]
for j in range(1,cJ-1):
    Flux += - k * (T[0,j] - T[1,j])/dy_C[j]
    Flux += - k * (T[-1,j] - T[-2,j])/dy_C[j]
    
print(Flux)

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

plt.figure()
plt.plot(nIter[1:-1],residuals[1:-1], "red")
plt.title('Residual convergence')
plt.xlabel('iterations')
plt.ylabel('residuals [-]')
plt.title('Residual')

plt.subplot(1,2,1)
plt.plot(yv[:,0], T[0,:])
plt.title('Temperature')
plt.xlabel('y [m]')
plt.ylabel('T [K]')
plt.show()


plt.subplot(1,2,2)
plt.plot(xv[0,:], T[:,0])
plt.title('Temperature')
plt.xlabel('x [m]')
plt.ylabel('T [K]')
plt.show()

## TASK 2 TODO
# 1. A convection term added => new equation 
# 2. New type of discretization: Hybrid scheme, Pe > 2 -> Upwind 
# 3. Use of either Gauss-Seidel or TDMA
# 4. Temperature in Kelvin!
# 5. Implement different termination criteria
#
#
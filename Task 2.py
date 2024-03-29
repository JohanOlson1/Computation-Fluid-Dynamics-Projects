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

def Interpolate_u(U_p, U_b, d_CV, d_N):
    f = 0.5 * d_CV / d_N
    U_interpolated = f * U_b + (1-f)*U_p     
    return U_interpolated

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
Cond_west   = 'Neumann'# Either 'Neumann' or 'Dirichlet'
caseID      =  15      # your case number to solve
k           =  1       # Conductivity
rho         =  1       # density
nIterations =  4000    # Max number of iterations
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
Flist = []

# Code
gamma = k/Cp                  # 
q_s = 100                     # Heat flow South Boundary   
Speed_factor = 1              # Factor to reduce the flow speed
T_New = 273.15 + 20           # Temperature for new Boundary Condition
 
## Diffusive and convective coefficient calculations
# Remark: No need to keep inside iteration lopp since independent of temp!
for i in range(1,cI-1):
	for j in range(1,cJ-1):
	    D[i,j,0] = gamma * dy_C[j] / dxe_C[i] # east diffusive
	    D[i,j,1] = gamma * dy_C[j] / dxw_C[i] # west diffusive
	    D[i,j,2] = gamma * dx_C[i] / dyn_C[j] # north diffusive
	    D[i,j,3] = gamma * dx_C[i] / dys_C[j] # south diffusive
	    
for i in range(2,cI-2):
	for j in range(2,cJ-2):
	    F[i,j,0] = rho * dy_C[j] * Interpolate_u(U[i,j], U[i+1,j], \
      dx_C[i], dxe_C[i])/Speed_factor    # east convective
	    F[i,j,1] = rho * dy_C[j] * Interpolate_u(U[i,j], U[i-1,j], \
      dx_C[i], dxw_C[i])/Speed_factor    # west convective
	    F[i,j,2] = rho * dx_C[i] * Interpolate_u(V[i,j], V[i,j+1], \
      dy_C[j], dyn_C[j])/Speed_factor    # north convective
	    F[i,j,3] = rho * dx_C[i] * Interpolate_u(V[i,j], V[i,j-1], \
      dy_C[j], dys_C[j])/Speed_factor    # south convective

# No need to interpolate velocity at boundaries!
for i in range(1,cI-1):
	F[i,-2,2] = rho * dx_C[i] * V[i,-1]/Speed_factor  # north convective
	F[i,1,3] = rho * dx_C[i] * V[i,0]/Speed_factor    # south convective
for j in range(1,cJ-1):
    F[-2,j,0] = rho * dy_C[j] * U[-1,j]/Speed_factor  # east convective
    F[1,j,1] = rho * dy_C[j] * U[0,j]/Speed_factor    # west convective
        
# Hybrid scheme coefficients calculations (taking into account boundary conditions afterwards!)
for i in range(1,cI-1):
	for j in range(1,cJ-1):
		coeffsT[i,j,0] = max([-F[i,j,0], D[i,j,0]-F[i,j,0]/2, 0])
		coeffsT[i,j,1] = max([F[i,j,1], D[i,j,1]+F[i,j,1]/2, 0])
		coeffsT[i,j,2] = max([-F[i,j,2], D[i,j,2]-F[i,j,2]/2, 0])
		coeffsT[i,j,3] = max([F[i,j,3], D[i,j,3]+F[i,j,3]/2, 0])
		coeffsT[i,j,4] = coeffsT[i,j,0] + coeffsT[i,j,1] + coeffsT[i,j,2] + \
        coeffsT[i,j,3] 

## East and West boundary
#Remark: Neumann Cond. -> coefficient is zero
for j in range(1, cJ-1):
    coeffsT[-2,j,0] = 0
    coeffsT[1,j,1] = 0
        
    # Recalculate a_p,
    coeffsT[-2,j,4] = coeffsT[-2,j,0] + coeffsT[-2,j,1] + coeffsT[-2,j,2] + \
    coeffsT[-2,j,3] 
    
    if Cond_west == 'Dirichlet':
        S_U[1,j] = max([F[1,j,1], D[1,j,1]+F[1,j,1]/2, 0]) * T_New 
        # Recalculate a_p
        coeffsT[1,j,4] = coeffsT[1,j,0] + coeffsT[1,j,1] + coeffsT[1,j,2] + \
    coeffsT[1,j,3] + max([F[1,j,1], D[1,j,1]+F[1,j,1]/2, 0])
    else:
        # Recalculate a_p
        coeffsT[1,j,4] = coeffsT[1,j,0] + coeffsT[1,j,1] + coeffsT[1,j,2] + \
        coeffsT[1,j,3]

# North and (South Boundary constant heat flow)
for i in range(1, cI-1):
    coeffsT[i,1,3] = 0
    coeffsT[i,-2,2] = 0
    
    S_U[i,1] = q_s * dx_C[i] /Cp
    S_U[i,-2] = max([-F[i,-2,2], D[i,-2,2]-F[i,-2,2]/2, 0]) * T_D
    
    # Recalculate a_p
    coeffsT[i,1,4] = coeffsT[i,1,0] + coeffsT[i,1,1] + coeffsT[i,1,2] + \
       coeffsT[i,1,3]
    
    # Recalculate a_p, actually not needed, but here for illustration! - S_P == a_N
    coeffsT[i,-2,4] = coeffsT[i,-2,0] + coeffsT[i,-2,1] + coeffsT[i,-2,2] + \
       coeffsT[i,-2,3] + max([-F[i,-2,2], D[i,-2,2]-F[i,-2,2]/2, 0])
    
# Temperature along west boundary
if Cond_west == 'Dirichlet':
    for j in range(cJ):
        T[0,j]  = T_New 

# Temperature along north boundary
for i in range(cI):
    T[i,-1] = T_D
    
for iter in range(nIterations): 
    
    # TDMA Solution
    if solver_type == 'TDMA':
        # Calculate all P, Q in for loop
        if (iter % 2) == 0:
            for j in range(1, cJ-1):
                for i in range(1, cI-1):
                    d_x[i] = coeffsT[i,j,2] * T[i,j+1] + coeffsT[i,j,3] \
                    * T[i,j-1] + S_U[i,j]
                    P_x[i] = coeffsT[i,j,0]/(coeffsT[i,j,4] - coeffsT[i,j,1]*P_x[i-1])
                    Q_x[i] = (d_x[i] + coeffsT[i,j,1] * Q_x[i-1])/\
                    (coeffsT[i,j,4] - coeffsT[i,j,1]*P_x[i-1])            
                # Backward differencing: East-West
                for k in range(cI-2, 0, -1):
                    T[k,j] = P_x[k] * T[k+1,j] + Q_x[k]
            
        else:      
            for i in range(1, cI-1):
                for j in range(1, cJ-1):
                    d_y[j] = coeffsT[i,j,0] * T[i+1,j] + coeffsT[i,j,1] * \
                    T[i-1,j] + S_U[i,j]
                    P_y[j] = coeffsT[i,j,2]/(coeffsT[i,j,4] - coeffsT[i,j,3]*P_y[j-1])
                    Q_y[j] = (d_y[j] + coeffsT[i,j,3] * Q_y[j-1])/(coeffsT[i,j,4] \
                       - coeffsT[i,j,3]*P_y[j-1])  
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
        T[i,0] = T[i,1]
    for j in range(cJ):
        if Cond_west == 'Neumann':
            T[0,j]  = T[1,j] 
        T[-1,j] = T[-2,j]
    
    # Compute residuals (taking into account normalization)
    r = 0
    for j in range(1,cJ-1):
        for i in range(1,cI-1):
            rhs = S_U[i,j] + coeffsT[i,j,0]*T[i+1,j] + coeffsT[i,j,1]*T[i-1,j] + \
            coeffsT[i,j,2]*T[i,j+1] + coeffsT[i,j,3]*T[i,j-1]
            r += np.abs(coeffsT[i,j,4]*T[i,j] - rhs)
            
    Inlet = 0; Outlet = 0
    for i in range(1,cI-1):
        # North Boundary
        Inlet += np.abs(rho * V[i,-1] * dx_C[i] * T[i,-1])
    for j in range(1,cJ-1):
        # East and West Boundary
        Outlet += np.abs(rho * U[0,j] * dy_C[j] * T[0,j]) + np.abs(rho * \
                        U[-1,j] * dy_C[j] * T[-1,j]) 
     
    r /= np.abs(Inlet - Outlet) # r/F
            
    Flist.append(np.abs(Inlet - Outlet))
    residuals.append(r) # fill with your residual value for the 
    nIter.append(iter)  # Append respective iter
    
    print('iteration: %d\nresT = %.5e\n\n' % (iter, residuals[-1]))
    
    # Check convergence
    if resTolerance>residuals[-1]:
        break
    
# Total Conservation
Inlet = Inlet*Cp
Outlet = Outlet*Cp
Flux_tot = 0
Diffusion_north = 0
Heat_flow_south = 0
for i in range(1,cI-1):
    Diffusion_north += - k * dx_C[i] * (T[i,-1] - T[i,-2])/dyn_C[-2]
    Heat_flow_south += q_s * dx_C[i] 

Flux_tot += Diffusion_north
Flux_tot += Heat_flow_south
Flux_tot += Inlet - Outlet
print("Netto flux in",Flux_tot)
print("South heat flow in", Heat_flow_south)
print("Diffusion in", Diffusion_north)
print("Flow in", Inlet)
print("Flow out", - Outlet)
print("Global Conservation",np.abs(Flux_tot/(Inlet + q_s + Diffusion_north)))

# Plotting 
xv, yv = np.meshgrid(xCoords_N, yCoords_N)

plt.figure()
plt.semilogy(nIter[1:-1],residuals[1:-1], "red", label="Residuals")
plt.semilogy(nIter[1:-1],Flist[1:-1], "blue", label="F")
plt.legend(loc="upper right")
plt.title('Residual convergence')
plt.xlabel('iterations')
plt.ylabel('Log(residuals)')
plt.title('Residuals and F')

plt.figure()
plt.plot(yv[1:-1,0], T[0,1:-1])
plt.title('Temperature Boundary West Boundary')
plt.xlabel('y [m]')
plt.ylabel('T [K]')
plt.show()

plt.figure()
mycmap1 = plt.get_cmap('coolwarm')
plt.contourf(xv, yv, T.T, cmap=mycmap1)
plt.colorbar()
plt.quiver(xv, yv, U.T, V.T, width=0.0011)
plt.title('Temperatures and Velocity Fields')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()

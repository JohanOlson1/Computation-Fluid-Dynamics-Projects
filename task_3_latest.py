# MTF072 Computational Fluid Dynamics
# Task 3: laminar lid-driven cavity
# Template prepared by:
# Gonzalo Montero Villar
# Department of Mechanics and Maritime Sciences
# Division of Fluid Dynamics
# villar@chalmers.se
# December 2020

def Rhie_Chow(Up, Ub, dA, a_p, ap_b, P_term):
    U_interpolated = 0.5*(Up + Ub) + (dA/(4*(a_p + ap_b)/2))*P_term
    return U_interpolated

#==============Packages needed=================
import matplotlib.pyplot as plt
import numpy as np

#================= Inputs =====================
### Inputs

# Fluid properties and B. C. inputs

UWall = 1 # velocity of the upper wall
rho   = 1 # density
nu    = 0.001 # kinematic viscosity

data_file = 'data_hybrid.txt' # data file where the given solution is stored

# Geometric inputs (fixed so that a fair comparison can be made)

mI = 11 # number of mesh points X direction. 
mJ = 11 # number of mesh points Y direction. 
xL =  1 # length in X direction
yL =  1 # length in Y direction

# Solver inputs

nIterations           =  25 # maximum number of iterations
n_inner_iterations_gs =  15 # amount of inner iterations when solving 
                              # pressure correction with Gauss-Seidel
resTolerance =  0.0001 # convergence criteria for residuals
                     # each variable
alphaUV =       0.5 # under relaxation factor for U and V
alphaP  =       0.5 # under relaxation factor for P

# ================ Code =======================
### Defintions

# For all the matrices the first input makes reference to the x coordinate
# and the second input to the y coordinate, (i+1) is east and (j+1) north

# Allocate all needed variables
nI = mI + 1                      # number of nodes in the X direction. nodes 
                                  # added in the boundaries
nJ = mJ + 1                      # number of nodes in the Y direction. nodes 
                                  # added in the boundaries
coeffsUV   = np.zeros((nI,nJ,5)) # coefficients for the U and V equation
                                  # E, W, N, S and P
sourceUV   = np.zeros((nI,nJ,2)) # source coefficients for the U and V equation
                                  # U and V
coeffsPp   = np.zeros((nI,nJ,5)) # coefficients for the pressure correction
                                  # equation E, W, N, S and P
sourcePp   = np.zeros((nI,nJ))   # source coefficients for the pressure
                                  # correction equation
U          = np.zeros((nI,nJ))   # U velocity matrix
V          = np.zeros((nI,nJ))   # V velocity matrix
P          = np.zeros((nI,nJ))   # pressure matrix
Pp         = np.zeros((nI,nJ))   # pressure correction matrix

massFlows  = np.zeros((nI,nJ,4)) # mass flows at the faces
                                  # m_e, m_w, m_n and m_s

residuals  = np.zeros((3,1))     # U, V and conitnuity residuals

# Generate mesh and compute geometric variables

# Allocate all variables matrices
xCoords_N = np.zeros((nI,nJ)) # X coords of the nodes
yCoords_N = np.zeros((nI,nJ)) # Y coords of the nodes
xCoords_M = np.zeros((mI,mJ)) # X coords of the mesh points
yCoords_M = np.zeros((mI,mJ)) # Y coords of the mesh points
dxe_N     = np.zeros((nI,nJ)) # X distance to east node
dxw_N     = np.zeros((nI,nJ)) # X distance to west node
dyn_N     = np.zeros((nI,nJ)) # Y distance to north node
dys_N     = np.zeros((nI,nJ)) # Y distance to south node
dx_CV      = np.zeros((nI,nJ)) # X size of the node
dy_CV      = np.zeros((nI,nJ)) # Y size of the node

residuals_U = []
residuals_V = []
residuals_c = []
nIter = []

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
            xCoords_N[i,j] = 0.5*(xCoords_M[i,j] + xCoords_M[i-1,j])
        if i == mI-1 and j>0:
            yCoords_N[i+1,j] = 0.5*(yCoords_M[i,j] + yCoords_M[i,j-1])
        if j > 0:
            yCoords_N[i,j] = 0.5*(yCoords_M[i,j] + yCoords_M[i,j-1])
        if j == mJ-1 and i>0:
            xCoords_N[i,j+1] = 0.5*(xCoords_M[i,j] + xCoords_M[i-1,j])

        # Fill dx_CV and dy_CV
        if i > 0:
            dx_CV[i,j] = xCoords_M[i,j] - xCoords_M[i-1,j]
        if j > 0:
            dy_CV[i,j] = yCoords_M[i,j] - yCoords_M[i,j-1]

xCoords_N[-1,:] = xL
yCoords_N[:,-1] = yL


# Fill dxe, dxw, dyn and dys
for i in range(1,nI-1):
    for j in range(1,nJ-1):
        dxe_N[i,j] = xCoords_N[i+1,j] - xCoords_N[i,j]
        dxw_N[i,j] = xCoords_N[i,j] - xCoords_N[i-1,j]
        dyn_N[i,j] = yCoords_N[i,j+1] - yCoords_N[i,j]
        dys_N[i,j] = yCoords_N[i,j] - yCoords_N[i,j-1]

# Initialize variable matrices

U[:,:] = 0
V[:,:] = 0
P[:,:] = 0

Dx = nu*rho/dx
Dy = nu*rho/dy

F_x = rho*dy
F_y = rho*dx

U[:,-1] = 1

# ================ Looping =======================
### Iteration start

for iter in range(nIterations):
    ### Boundary Conditions
    ## Impose boundary conditions for velocities, only the top boundary wall
        # TO-DO: Make it implicit instead of explicit?
    #U[:,-1] = 1
    
    # Impose pressure boundary condition, all homogeneous Neumann
    # How about no
    
    ### U,V Coefficents
    # First, north and south boundaries
    ### Question: U and V in all 
    for i in range(2,nI-2):
        # North
        j = -2
        coeffsUV[i,j,0] = Dx + max([0,- F_x*(U[i,j]+U[i+1,j])/2])
        coeffsUV[i,j,1] = Dx + max([F_x*(U[i,j]+U[i-1,j])/2,0])
        coeffsUV[i,j,2] = 0
        coeffsUV[i,j,3] = Dy + max([F_y*(V[i,j]+V[i,j-1])/2,0])
        coeffsUV[i,j,4] = coeffsUV[i,j,0] + coeffsUV[i,j,1] + coeffsUV[i,j,2] + coeffsUV[i,j,3] + (2*Dy + max([F_y*(V[i,j+1]),0]))
        sourceUV[i,j,0] = (P[i-1,j]-P[i+1,j])*dy + (2*Dy + max([F_y*(V[i,j+1]),0]))*U[i,j+1]
        sourceUV[i,j,1] = ((P[i,j-1]+P[i,j])/2-P[i,j+1])*dx
        
        # South
        j = 1
        coeffsUV[i,j,0] = Dx + max([0,- F_x*(U[i,j]+U[i+1,j])/2])
        coeffsUV[i,j,1] = Dx + max([F_x*(U[i,j]+U[i-1,j])/2,0])
        coeffsUV[i,j,2] = Dy + max([0,- F_y*(V[i,j]+V[i,j+1])/2])
        coeffsUV[i,j,3] = 0
        coeffsUV[i,j,4] = coeffsUV[i,j,0] + coeffsUV[i,j,1] + coeffsUV[i,j,2] + coeffsUV[i,j,3]  
        sourceUV[i,j,0] = (P[i-1,j]-P[i+1,j])*dy/2
        sourceUV[i,j,1] = (P[i,j-1] - (P[i,j] + P[i,j+1])/2)*dx

    # Second, east and west boundaries
    for j in range(2,nJ-2):
        # East
        i = -2
        coeffsUV[i,j,0] = 0
        coeffsUV[i,j,1] = Dx + max([F_x*(U[i,j]+U[i-1,j])/2,0])
        coeffsUV[i,j,2] = Dy + max([0,- F_y*(V[i,j]+V[i,j+1])/2])
        coeffsUV[i,j,3] = Dy + max([F_y*(V[i,j]+V[i,j-1])/2,0])
        coeffsUV[i,j,4] = coeffsUV[i,j,0] + coeffsUV[i,j,1] + coeffsUV[i,j,2] + coeffsUV[i,j,3]  
        sourceUV[i,j,0] = ((P[i-1,j]+P[i,j])/2-P[i+1,j])*dy
        sourceUV[i,j,1] = (P[i,j-1]-P[i,j+1])*dx/2
        
        # West
        i = 1
        coeffsUV[i,j,0] = Dx + max([0,- F_x*(U[i,j]+U[i+1,j])/2])
        coeffsUV[i,j,1] = 0
        coeffsUV[i,j,2] = Dy + max([0,- F_y*(V[i,j]+V[i,j+1])/2])
        coeffsUV[i,j,3] = Dy + max([F_y*(V[i,j]+V[i,j-1])/2,0])
        coeffsUV[i,j,4] = coeffsUV[i,j,0] + coeffsUV[i,j,1] + coeffsUV[i,j,2] + coeffsUV[i,j,3] 
        sourceUV[i,j,0] = (P[i-1,j]-(P[i+1,j]+P[i,j]/2))*dy
        sourceUV[i,j,1] = (P[i,j-1]-P[i,j+1])*dx/2
         
    ## Compute coefficients at corner nodes (one step inside)
    
    # North-West Corner
    i = 1; j = -2;
    coeffsUV[i,j,0] = Dx + max([0,- F_x*(U[i,j]+U[i+1,j])/2])
    coeffsUV[i,j,1] = 0
    coeffsUV[i,j,2] = 0
    coeffsUV[i,j,3] = Dy + max([F_y*(V[i,j]+V[i,j-1])/2,0])
    coeffsUV[i,j,4] = coeffsUV[i,j,0] + coeffsUV[i,j,1] + coeffsUV[i,j,2] + coeffsUV[i,j,3] + (2*Dy + max([F_y*(V[i,j+1]),0]))
    sourceUV[i,j,0] = (P[i-1,j]-(P[i+1,j]+P[i,j]/2))*dy + (2*Dy + max([F_y*(V[i,j+1]),0]))*U[i,j+1]
    sourceUV[i,j,1] = ((P[i,j-1]+P[i,j])/2-P[i,j+1])*dx 
    
    # North-East Corner
    i = -2; j = -2;
    coeffsUV[i,j,0] = 0
    coeffsUV[i,j,1] = Dx + max([F_x*(U[i,j]+U[i-1,j])/2,0])
    coeffsUV[i,j,2] = 0
    coeffsUV[i,j,3] = Dy + max([F_y*(V[i,j]+V[i,j-1])/2,0])
    coeffsUV[i,j,4] = coeffsUV[i,j,0] + coeffsUV[i,j,1] + coeffsUV[i,j,2] + coeffsUV[i,j,3] + (2*Dy + max([F_y*(V[i,j+1]),0])) 
    sourceUV[i,j,0] = ((P[i-1,j]+P[i,j])/2-P[i+1,j])*dy + (2*Dy + max([F_y*(V[i,j+1]),0]))*U[i,j+1]
    sourceUV[i,j,1] = ((P[i,j-1]+P[i,j])/2-P[i,j+1])*dx
    
    # South-West Corner
    i = 1; j = 1; 
    coeffsUV[i,j,0] = Dx + max([0,- F_x*(U[i,j]+U[i+1,j])/2])
    coeffsUV[i,j,1] = 0
    coeffsUV[i,j,2] = Dy + max([0,- F_y*(V[i,j]+V[i,j+1])/2])
    coeffsUV[i,j,3] = 0
    coeffsUV[i,j,4] = coeffsUV[i,j,0] + coeffsUV[i,j,1] + coeffsUV[i,j,2] + coeffsUV[i,j,3] 
    sourceUV[i,j,0] = (P[i-1,j]-(P[i+1,j]+P[i,j]/2))*dy
    sourceUV[i,j,1] = (P[i,j-1] - (P[i,j] + P[i,j+1])/2)*dx
    
    # South-East Corner
    i = -2; j = 1;
    coeffsUV[i,j,0] = Dx + max([0,- F_x*(U[i,j]+U[i+1,j])/2])
    coeffsUV[i,j,1] = 0
    coeffsUV[i,j,2] = Dy + max([0,- F_y*(V[i,j]+V[i,j+1])/2])
    coeffsUV[i,j,3] = 0
    coeffsUV[i,j,4] = coeffsUV[i,j,0] + coeffsUV[i,j,1] + coeffsUV[i,j,2] + coeffsUV[i,j,3] 
    sourceUV[i,j,0] = ((P[i-1,j]+P[i,j])/2-P[i+1,j])*dy
    sourceUV[i,j,1] = (P[i,j-1] - (P[i,j] + P[i,j+1])/2)*dx
    
    
    ## Compute coefficients for inner nodes
    for i in range(2,nI-2):
        for j in range(2,nJ-2):
            coeffsUV[i,j,0] = Dx + max([0,- F_x*(U[i,j]+U[i+1,j])/2])
            coeffsUV[i,j,1] = Dx + max([F_x*(U[i,j]+U[i-1,j])/2,0])
            coeffsUV[i,j,2] = Dy + max([0,- F_y*(V[i,j]+V[i,j+1])/2])
            coeffsUV[i,j,3] = Dy + max([F_y*(V[i,j]+V[i,j-1])/2,0])
            coeffsUV[i,j,4] = coeffsUV[i,j,0] + coeffsUV[i,j,1] + coeffsUV[i,j,2] + coeffsUV[i,j,3]  
            sourceUV[i,j,0] = (P[i-1,j]-P[i+1,j])*dy/2
            sourceUV[i,j,1] = (P[i,j-1]-P[i,j+1])*dx/2
        
    ### Implicit Under-Relaxation 
    ### SIMPLE - 1. U*,V*
    ## Solve for U* and V* using Gauss-Seidel
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            U[i,j] = alphaUV*(coeffsUV[i,j,0] * U[i+1,j] + coeffsUV[i,j,1] * U[i-1,j] \
                    + coeffsUV[i,j,2] * U[i,j+1] + coeffsUV[i,j,3] * U[i,j-1] + sourceUV[i,j,0])/coeffsUV[i,j,4] \
                + (1-alphaUV)*U[i,j]
            V[i,j] = alphaUV*(coeffsUV[i,j,0] * V[i+1,j] + coeffsUV[i,j,1] * V[i-1,j] \
                    + coeffsUV[i,j,2] * V[i,j+1] + coeffsUV[i,j,3] * V[i,j-1] + sourceUV[i,j,1])/coeffsUV[i,j,4] \
                + (1-alphaUV)*V[i,j]
                  
    ## Calculate at the faces using Rhie-Chow for the face velocities
    # First, North and South
    for i in range(2,nI-2):
        # North
        j = -2;
        massFlows[i,j,0] = Rhie_Chow(U[i,j], U[i+1,j], dy, coeffsUV[i,j,4], \
            coeffsUV[i+1,j,4], P[i+2,j] - 3*P[i+1,j] + 3*P[i,j] + P[i-1,j])
        massFlows[i,j,1] = Rhie_Chow(U[i,j], U[i-1,j], dy, coeffsUV[i,j,4], \
            coeffsUV[i-1,j,4], P[i+1,j] - 3*P[i,j] + 3*P[i-1,j] + P[i-2,j])
        massFlows[i,j,2] = V[i,j+1]
        massFlows[i,j,3] = Rhie_Chow(V[i,j], V[i,j-1], dx, coeffsUV[i,j,4], \
            coeffsUV[i,j-1,4], P[i,j+1] - 3*P[i,j] + 3*P[i,j-1] + P[i,j-2]) 
                
        # South
        j = 1;
        massFlows[i,j,0] = Rhie_Chow(U[i,j], U[i+1,j], dy, coeffsUV[i,j,4], \
            coeffsUV[i+1,j,4], P[i+2,j] - 3*P[i+1,j] + 3*P[i,j] + P[i-1,j])
        massFlows[i,j,1] = Rhie_Chow(U[i,j], U[i-1,j], dy, coeffsUV[i,j,4], \
            coeffsUV[i-1,j,4], P[i+1,j] - 3*P[i,j] + 3*P[i-1,j] + P[i-2,j])
        massFlows[i,j,2] = Rhie_Chow(V[i,j], V[i,j+1], dx, coeffsUV[i,j,4], \
            coeffsUV[i,j+1,4], P[i,j+2] - 3*P[i,j+1] + 3*P[i,j] + P[i,j-1])
        massFlows[i,j,3] = V[i,j-1] 
    
    # Second, East and West
    for j in range(2,nI-2):
        # East
        i = -2;
        massFlows[i,j,0] = U[i+1,j]
        massFlows[i,j,1] = Rhie_Chow(U[i,j], U[i-1,j], dy, coeffsUV[i,j,4], \
            coeffsUV[i-1,j,4], P[i+1,j] - 3*P[i,j] + 3*P[i-1,j] + P[i-2,j])
        massFlows[i,j,2] = Rhie_Chow(V[i,j], V[i,j+1], dx, coeffsUV[i,j,4], \
            coeffsUV[i,j+1,4], P[i,j+2] - 3*P[i,j+1] + 3*P[i,j] + P[i,j-1])
        massFlows[i,j,3] = Rhie_Chow(V[i,j], V[i,j-1], dx, coeffsUV[i,j,4], \
            coeffsUV[i,j-1,4], P[i,j+1] - 3*P[i,j] + 3*P[i,j-1] + P[i,j-2]) 
            
        # West
        i = 1;
        massFlows[i,j,0] = Rhie_Chow(U[i,j], U[i+1,j], dy, coeffsUV[i,j,4], \
            coeffsUV[i+1,j,4], P[i+2,j] - 3*P[i+1,j] + 3*P[i,j] + P[i-1,j])
        massFlows[i,j,1] = U[i-1,j]
        massFlows[i,j,2] = Rhie_Chow(V[i,j], V[i,j+1], dx, coeffsUV[i,j,4], \
            coeffsUV[i,j+1,4], P[i,j+2] - 3*P[i,j+1] + 3*P[i,j] + P[i,j-1])
        massFlows[i,j,3] = Rhie_Chow(V[i,j], V[i,j-1], dx, coeffsUV[i,j,4], \
            coeffsUV[i,j-1,4], P[i,j+1] - 3*P[i,j] + 3*P[i,j-1] + P[i,j-2]) 
    
    # North-West Corner
    i = 1; j = -2;
    massFlows[i,j,0] = Rhie_Chow(U[i,j], U[i+1,j], dy, coeffsUV[i,j,4], \
        coeffsUV[i+1,j,4], P[i+2,j] - 3*P[i+1,j] + 3*P[i,j] + P[i-1,j])
    massFlows[i,j,1] = U[i-1,j]
    massFlows[i,j,2] = V[i,j+1]
    massFlows[i,j,3] = Rhie_Chow(V[i,j], V[i,j-1], dx, coeffsUV[i,j,4], \
        coeffsUV[i,j-1,4], P[i,j+1] - 3*P[i,j] + 3*P[i,j-1] + P[i,j-2]) 
    
    # North-East Corner
    i = -2; j = -2; 
    massFlows[i,j,0] = U[i+1,j]
    massFlows[i,j,1] = Rhie_Chow(U[i,j], U[i-1,j], dy, coeffsUV[i,j,4], \
        coeffsUV[i-1,j,4], P[i+1,j] - 3*P[i,j] + 3*P[i-1,j] + P[i-2,j])
    massFlows[i,j,2] = V[i,j+1]
    massFlows[i,j,3] = Rhie_Chow(V[i,j], V[i,j-1], dx, coeffsUV[i,j,4], \
        coeffsUV[i,j-1,4], P[i,j+1] - 3*P[i,j] + 3*P[i,j-1] + P[i,j-2]) 
    
    # South-West Corner
    i = 1; j = 1;
    massFlows[i,j,0] = Rhie_Chow(U[i,j], U[i+1,j], dy, coeffsUV[i,j,4], \
        coeffsUV[i+1,j,4], P[i+2,j] - 3*P[i+1,j] + 3*P[i,j] + P[i-1,j])
    massFlows[i,j,1] = U[i-1,j]
    massFlows[i,j,2] = Rhie_Chow(V[i,j], V[i,j+1], dx, coeffsUV[i,j,4], \
        coeffsUV[i,j+1,4], P[i,j+2] - 3*P[i,j+1] + 3*P[i,j] + P[i,j-1])
    massFlows[i,j,3] = V[i,j-1]
    
    # South-East Corner
    i = -2; j = 1;
    massFlows[i,j,0] = U[i+1,j]
    massFlows[i,j,1] = Rhie_Chow(U[i,j], U[i-1,j], dy, coeffsUV[i,j,4], \
        coeffsUV[i-1,j,4], P[i+1,j] - 3*P[i,j] + 3*P[i-1,j] + P[i-2,j])
    massFlows[i,j,2] = Rhie_Chow(V[i,j], V[i,j+1], dx, coeffsUV[i,j,4], \
        coeffsUV[i,j+1,4], P[i,j+2] - 3*P[i,j+1] + 3*P[i,j] + P[i,j-1])
    massFlows[i,j,3] = V[i,j-1] 
    
    # Inner Nodes
    for i in range(2,nI-2):
        for j in range(2,nJ-2): 
            massFlows[i,j,0] = Rhie_Chow(U[i,j], U[i+1,j], dy, coeffsUV[i,j,4], \
                coeffsUV[i+1,j,4], P[i+2,j] - 3*P[i+1,j] + 3*P[i,j] + P[i-1,j])
            massFlows[i,j,1] = Rhie_Chow(U[i,j], U[i-1,j], dy, coeffsUV[i,j,4], \
                coeffsUV[i-1,j,4], P[i+1,j] - 3*P[i,j] + 3*P[i-1,j] + P[i-2,j])
            massFlows[i,j,2] = Rhie_Chow(V[i,j], V[i,j+1], dx, coeffsUV[i,j,4], \
                coeffsUV[i,j+1,4], P[i,j+2] - 3*P[i,j+1] + 3*P[i,j] + P[i,j-1])
            massFlows[i,j,3] = Rhie_Chow(V[i,j], V[i,j-1], dx, coeffsUV[i,j,4], \
                coeffsUV[i,j-1,4], P[i,j+1] - 3*P[i,j] + 3*P[i,j-1] + P[i,j-2]) 
    
    ### SIMPLE - 2. Pressure Correction
    ## Calculate pressure correction equation coefficients 
    Pp_x = rho*dy*dy*alphaUV
    Pp_y = rho*dx*dx*alphaUV
    # First, North and South
    for i in range(2,nI-2):
        # North
        j = -2
        coeffsPp[i,j,0] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i+1,j,4])/2)
        coeffsPp[i,j,1] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i-1,j,4])/2)
        coeffsPp[i,j,2] = 0
        coeffsPp[i,j,3] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j-1,4])/2)
        coeffsPp[i,j,4] = coeffsPp[i,j,0] + coeffsPp[i,j,1] + coeffsPp[i,j,2] + coeffsPp[i,j,3]
        sourcePp[i,j]   = - rho*(massFlows[i,j,0] - massFlows[i,j,1])*dy - \
            rho*(massFlows[i,j,2] - massFlows[i,j,3])*dx 
        
        # South
        j = 1
        coeffsPp[i,j,0] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i+1,j,4])/2)
        coeffsPp[i,j,1] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i-1,j,4])/2)
        coeffsPp[i,j,2] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j+1,4])/2)
        coeffsPp[i,j,3] = 0
        coeffsPp[i,j,4] = coeffsPp[i,j,0] + coeffsPp[i,j,1] + coeffsPp[i,j,2] + coeffsPp[i,j,3]
        sourcePp[i,j]   = - rho*(massFlows[i,j,0] - massFlows[i,j,1])*dy - \
            rho*(massFlows[i,j,2] - massFlows[i,j,3])*dx 
        
    # Second, East and West
    for j in range(2,nJ-2):
        # East
        i = -2
        coeffsPp[i,j,0] = 0
        coeffsPp[i,j,1] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i-1,j,4])/2)
        coeffsPp[i,j,2] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j+1,4])/2)
        coeffsPp[i,j,3] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j-1,4])/2)
        coeffsPp[i,j,4] = coeffsPp[i,j,0] + coeffsPp[i,j,1] + coeffsPp[i,j,2] + coeffsPp[i,j,3]
        sourcePp[i,j]   = - rho*(massFlows[i,j,0] - massFlows[i,j,1])*dy - \
            rho*(massFlows[i,j,2] - massFlows[i,j,3])*dx 
            
        # West
        i = 1
        coeffsPp[i,j,0] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i+1,j,4])/2)
        coeffsPp[i,j,1] = 0
        coeffsPp[i,j,2] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j+1,4])/2)
        coeffsPp[i,j,3] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j-1,4])/2)
        coeffsPp[i,j,4] = coeffsPp[i,j,0] + coeffsPp[i,j,1] + coeffsPp[i,j,2] + coeffsPp[i,j,3]
        sourcePp[i,j]   = - rho*(massFlows[i,j,0] - massFlows[i,j,1])*dy - \
            rho*(massFlows[i,j,2] - massFlows[i,j,3])*dx 
        
    # North-West Corner
    i = 1; j = -2;
    coeffsPp[i,j,0] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i+1,j,4])/2)
    coeffsPp[i,j,1] = 0
    coeffsPp[i,j,2] = 0
    coeffsPp[i,j,3] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j-1,4])/2)
    coeffsPp[i,j,4] = coeffsPp[i,j,0] + coeffsPp[i,j,1] + coeffsPp[i,j,2] + coeffsPp[i,j,3]
    sourcePp[i,j]   = - rho*(massFlows[i,j,0] - massFlows[i,j,1])*dy - \
        rho*(massFlows[i,j,2] - massFlows[i,j,3])*dx 
        
    # North-East Corner
    i = -2; j = -2;
    coeffsPp[i,j,0] = 0
    coeffsPp[i,j,1] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i-1,j,4])/2)
    coeffsPp[i,j,2] = 0
    coeffsPp[i,j,3] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j-1,4])/2)
    coeffsPp[i,j,4] = coeffsPp[i,j,0] + coeffsPp[i,j,1] + coeffsPp[i,j,2] + coeffsPp[i,j,3]
    sourcePp[i,j]   = - rho*(massFlows[i,j,0] - massFlows[i,j,1])*dy - \
        rho*(massFlows[i,j,2] - massFlows[i,j,3])*dx 
        
    # South-West
    i = 1; j = 1;
    coeffsPp[i,j,0] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i+1,j,4])/2)
    coeffsPp[i,j,1] = 0
    coeffsPp[i,j,2] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j+1,4])/2)
    coeffsPp[i,j,3] = 0
    coeffsPp[i,j,4] = coeffsPp[i,j,0] + coeffsPp[i,j,1] + coeffsPp[i,j,2] + coeffsPp[i,j,3]
    sourcePp[i,j]   = - rho*(massFlows[i,j,0] - massFlows[i,j,1])*dy - \
        rho*(massFlows[i,j,2] - massFlows[i,j,3])*dx 
    
    # South-East
    i = -2; j = 1;
    coeffsPp[i,j,0] = 0
    coeffsPp[i,j,1] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i-1,j,4])/2)
    coeffsPp[i,j,2] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j+1,4])/2)
    coeffsPp[i,j,3] = 0
    coeffsPp[i,j,4] = coeffsPp[i,j,0] + coeffsPp[i,j,1] + coeffsPp[i,j,2] + coeffsPp[i,j,3]
    sourcePp[i,j]   = - rho*(massFlows[i,j,0] - massFlows[i,j,1])*dy - \
        rho*(massFlows[i,j,2] - massFlows[i,j,3])*dx 
    
    # Inner nodes
    for i in range(2,nI-2):
        for j in range(2,nJ-2):
            coeffsPp[i,j,0] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i+1,j,4])/2)
            coeffsPp[i,j,1] = Pp_x/((coeffsUV[i,j,4]+coeffsUV[i-1,j,4])/2)
            coeffsPp[i,j,2] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j+1,4])/2)
            coeffsPp[i,j,3] = Pp_y/((coeffsUV[i,j,4]+coeffsUV[i,j-1,4])/2)
            coeffsPp[i,j,4] = coeffsPp[i,j,0] + coeffsPp[i,j,1] + coeffsPp[i,j,2] + coeffsPp[i,j,3]
            sourcePp[i,j]   = - (massFlows[i,j,0] - massFlows[i,j,1])*rho*dy - \
                (massFlows[i,j,2] - massFlows[i,j,3])*rho*dx 
            
       
    # Solve for pressure correction (Note that more that one loop is used)
    for iter_gs in range(n_inner_iterations_gs):
        for j in range(1,nJ-1):
            for i in range(1,nI-1):    
                Pp[i,j] = (coeffsPp[i,j,0] * Pp[i+1,j] + coeffsPp[i,j,1] * Pp[i-1,j] \
                      + coeffsPp[i,j,2] * Pp[i,j+1] + coeffsPp[i,j,3] * Pp[i,j-1] + sourcePp[i,j])/coeffsPp[i,j,4]
               
    # Set Pp with reference to node (2,2) and copy to boundaries
    Pp[1:-1,1:-1] -= Pp[1,1]
    
    ### SIMPLE - 3. Correct U,V,P
    # Correct velocities, pressure and mass flows, : Mass flow needed?
    for i in range(1,nI-1):
        for j in range(1,nJ-1): 
            P[i,j] += alphaP*Pp[i,j] 
            U[i,j] += (dy*alphaUV*(Pp[i-1,j] - Pp[i+1,j])/2)/((coeffsUV[i,j,4] + coeffsUV[i,j,4])/2)
            V[i,j] += (dx*alphaUV*(Pp[i,j-1] - Pp[i,j+1])/2)/((coeffsUV[i,j,4] + coeffsUV[i,j,4])/2)
            
            ### Correct mass flow here?
            # (massFlows[5,5,0] - massFlows[5,5,1]) + (massFlows[5,5,2] - massFlows[5,5,3])
    
    # impose zero mass flow at the boundaries :  
    V[:,-1] = 0 # North
    U[:,0]  = 0 # U South
    V[:,0]  = 0 # V South
    U[-1,:] = 0 # U East
    V[-1,:] = 0 # V East
    U[0,:]  = 0 # U West
    V[0,:]  = 0 # U West
        
    # Copy P to boundaries
    P[0,:]  = P[1,:]
    P[-1,:] = P[-2,:]
    P[:,0]  = P[:,1]
    P[:,-1] = P[:,-2]
    
    ### Compute Residuals
    
    rU = 0
    rV = 0
    for j in range(1,nJ-1):
        for i in range(1,nI-1):
            rhsU = coeffsUV[i,j,0]*U[i-1,j] + coeffsUV[i,j,1]*U[i+1,j] + \
            coeffsUV[i,j,2]*U[i,j-1] + coeffsUV[i,j,3]*U[i,j+1]
            rU += np.abs(coeffsUV[i,j,4]*U[i,j] - rhsU)
            
            rhsV = coeffsUV[i,j,0]*V[i-1,j] + coeffsUV[i,j,1]*V[i+1,j] + \
            coeffsUV[i,j,2]*V[i,j-1] + coeffsUV[i,j,3]*V[i,j+1]
            rV += np.abs(coeffsUV[i,j,4]*V[i,j] - rhsV)
            
    residuals_U.append(rU) # U momentum residual
    residuals_V.append(rV) # V momentum residual
    residuals_c.append(1) # continuity residual

# =============================================================================
#     residuals_U[-1] = rU
#     residuals_V[-1] = rV
#     residuals_c[-1] = 1
# =============================================================================
    nIter.append(iter)

    print('iteration: %d\nresU = %.5e, resV = %.5e, resCon = %.5e\n\n'\
        % (iter, residuals_U[-1], residuals_V[-1], residuals_c[-1]))
    
    #  Check convergence
    if resTolerance>max([residuals_U[-1], residuals_V[-1], residuals_c[-1]]):
        break

### Plotting
xv, yv = np.meshgrid(xCoords_N, yCoords_N)
# Plot mesh
plt.figure()
plt.plot(xCoords_N, yCoords_N, ".", markersize=1, color = "black")
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Computational mesh')

# Plot results

plt.figure()

# U velocity contour
ax = plt.subplot(2,3,1)
plt.title('U velocity [m/s]')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
mycmap1 = plt.get_cmap('coolwarm')
cp = ax.contourf(xCoords_N, yCoords_N, U, cmap=mycmap1)
plt.colorbar(cp)

# V velocity contour
ax2 = plt.subplot(2,3,2)
plt.title('V velocity [m/s]')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
cp2 = ax2.contourf(xCoords_N, yCoords_N, V, cmap=mycmap1)
plt.colorbar(cp2)

# P contour
ax3 = plt.subplot(2,3,3)
plt.title('Pressure [Pa]')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
cp3 = ax3.contourf(xCoords_N, yCoords_N, P, cmap=mycmap1)
plt.colorbar(cp3)

# Vector plot
plt.subplot(2,3,4)
plt.title('Vector plot of the velocity field')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.quiver(xCoords_N,yCoords_N,U,V)

# Residuals
plt.subplot(2,3,6)
plt.plot(nIter[1:-1],residuals_U[1:-1], "red")
plt.plot(nIter[1:-1],residuals_V[1:-1], "green")
plt.plot(nIter[1:-1],residuals_c[1:-1], "yellow")
plt.title('Residual convergence')
plt.xlabel('iterations')
plt.ylabel('residuals [-]')
#plt.legend('U momentum','V momentum', 'Continuity')
plt.title('Residuals')
plt.show()

# =============================================================================
# # Comparison with data
# data=np.genfromtxt(data_file, skip_header=1)
# uInterp = np.zeros((nJ-2,1))
# vInterp = np.zeros((nJ-2,1))
# for j in range(1,nJ-1):
#     for i in range(1,nI-1):
#         if xCoords_N[i,j]<0.5 and xCoords_N[i+1,j]>0.5:
#             uInterp[j-1] = (U[i+1,j] + U[i,j])*0.5
#             vInterp[j-1] = (V[i+1,j] + V[i,j])*0.5
#             break
#         elif abs(xCoords_N[i,j]-0.5) < 0.000001:
#             uInterp[j-1] = U[i,j]
#             vInterp[j-1] = V[i,j]
#             break
# 
# plt.subplot(2,3,5)
# plt.plot(data[:,0],data[:,2],'r.',markersize=20,label='data U')
# plt.plot(data[:,1],data[:,2],'b.',markersize=20,label='data V')
# plt.plot(uInterp,yCoords_N[1,1:-1],'k',label='sol U')
# plt.plot(vInterp,yCoords_N[1,1:-1],'g',label='sol V')
# plt.title('Comparison with data at x = 0.5')
# plt.xlabel('u, v [m/s]')
# plt.ylabel('y [m]')
# plt.legend()
# 
# plt.subplot(2,3,6)
# plt.title('Residual convergence')
# plt.xlabel('iterations')
# plt.ylabel('residuals [-]')
# #plt.legend('U momentum','V momentum', 'Continuity')
# plt.title('Residuals')
# plt.show()
# =============================================================================



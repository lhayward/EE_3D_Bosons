import numpy as np
#from scipy import linalg

#Ev_tol = 10e-16

#########################################  omega_3d  ##########################################
# Calculates omega for given kx, ky, kz and mass 
###############################################################################################
def omega_3d(kx, ky, kz, mass):
  return np.sqrt( 4*np.sin(kx/2.0)**2 + 4*np.sin(ky/2.0)**2 + 4*np.sin(kz/2.0)**2 + mass**2 )

#####################################  getCorrelators_3d  #####################################
# Calculates the phi-phi and pi-pi correlator for a 3D free bosonic system
# Input paramters: 
#    L: length of lattice (integer)
#    bc_x: boundary conditions along x ('PBC' or 'APBC')
#    bc_y: boundary conditions along y ('PBC' or 'APBC')
#    bc_z: boundary conditions along y ('PBC' or 'APBC')
#    mass: boson mass (float)
#    r, rprime: (x,y,z) coordinates on lattice of the two phi/pi variables (floats)
###############################################################################################
def getCorrelators_3d(L, bc_x, bc_y, bc_z, mass, r, rprime):
  d = 3
  bc_err = False
  (x,y,z) = r
  (xp, yp, zp) = rprime
  
  #BC along x:
  if bc_x == 'PBC':
    kx = (2*np.array(range(0,L)))*np.pi/L
  elif bc_x == 'APBC':
    kx = ( 2*np.array(range(0,L)) + 1)*np.pi/L
  else:
    kx = np.zeros(L)
    print "*** Boundary condition %s along x is not supported ***" %bc_x
    bc_err = True 
  
  #BC along y:
  if bc_y == 'PBC':
    ky = (2*np.array(range(0,L)))*np.pi/L
  elif bc_y == 'APBC':
    ky = ( 2*np.array(range(0,L)) + 1)*np.pi/L 
  else:
    kx = np.zeros(L)
    print "*** Boundary condition %s along y is not supported ***" %bc_y
    bc_err = True 
  
  #BC along z:
  if bc_z == 'PBC':
    kz = (2*np.array(range(0,L)))*np.pi/L
  elif bc_z == 'APBC':
    kz = ( 2*np.array(range(0,L)) + 1)*np.pi/L 
  else:
    kz = np.zeros(L)
    print "*** Boundary condition %s along z is not supported ***" %bc_z
    bc_err = True 
  
  phiphi = 0
  pipi   = 0
  if not bc_err:
    for kyy in ky:
      for kzz in kz:
        omega = omega_3d(kx,kyy,kzz,mass)
        phiphi = phiphi + sum( np.cos(kx*(x-xp))*np.cos(kyy*(y-yp))*np.cos(kzz*(z-zp))/omega )
        pipi   = pipi   + sum( np.cos(kx*(x-xp))*np.cos(kyy*(y-yp))*np.cos(kzz*(z-zp))*omega )

  return phiphi/(2*L**d), pipi/(2*L**d)

#####################################  getXAPA_edge_corr  #####################################
###############################################################################################
def getXAPA_edge_corr(L, LA_x, LA_y, LA_z, bc_x, bc_y, bc_z, mass):

  #Calculate all needing correlators:
  X_from0 = np.zeros((LA_x,LA_y,LA_z))
  P_from0 = np.zeros((LA_x,LA_y,LA_z))
  for x in range(LA_x):
    for y in range(LA_y):
      for z in range(LA_z):
        X_from0[x,y,z], P_from0[x,y,z] = getCorrelators_3d(L, bc_x, bc_y, bc_z, mass, (0,0,0), (x,y,z))
  
  sitesA = np.zeros(LA_x*LA_y*LA_z).tolist() #the indices of the sites in region A
  count = 0
  for x in range(0,LA_x):
    for y in range(0,LA_y):
      for z in range(0,LA_z):
        sitesA[count] = site((L,L,L),x,y,z)
        count = count + 1
  
  #Calculate XA and PA:
  NA = len(sitesA)
  XA = np.zeros((NA,NA))
  PA = np.zeros((NA,NA))
  for iA, sitei in enumerate(sitesA):
    for jA, sitej in enumerate(sitesA): 
      xi = sitei%L #sitei%Lx
      yi = ((sitei-xi)/L)%L #((sitei-xi)/Lx)%Ly
      zi = (sitei-xi-yi*L)/(L*L) #(sitei-x-y*Lx)/(Lx*Ly)
    
      xj = sitej%L 
      yj = ((sitej-xj)/L)%L 
      zj = (sitej-xj-yj*L)/(L*L) 
    
      XA[iA,jA] = X_from0[abs(xj-xi),abs(yj-yi),abs(zj-zi)]
      PA[iA,jA] = P_from0[abs(xj-xi),abs(yj-yi),abs(zj-zi)]
  
  return XA, PA

######################################  getXAPA_edge_K  #######################################
###############################################################################################
def getXAPA_edge_K(L, LA_x, LA_y, LA_z, bc_x, bc_y, bc_z, mass):
  d = 3
  Lx=Ly=Lz=L
  Ltup = (Lx,Ly,Lz)
  Ns = Lx*Ly*Lz
  K = np.zeros((Ns, Ns))
  
  XA = 0
  PA = 0
  
  bc_list = ['PBC', 'APBC', 'OBCNeu', 'OBCDir']
  
  if (bc_x in bc_list) and (bc_y in bc_list) and (bc_z in bc_list):
    #Loop to calculate the matrix K:
    for x in range(Lx):
      for y in range(Ly):
        for z in range(Lz):
          site_xyz = site(Ltup,x,y,z)
          xp = (x+1)%Lx
          yp = (y+1)%Ly
          zp = (z+1)%Lz
        
          #Onsite term:
          K[site_xyz,site_xyz] = 6.0 + (mass ** (2))
          #Adjust onsite term for Neumann BC (on boundaries):
          if bc_x == 'OBCNeu' and (xp < x or x==0):
            K[site_xyz,site_xyz] = K[site_xyz,site_xyz] - 1
          if bc_y == 'OBCNeu' and (yp < y or y==0):
            K[site_xyz,site_xyz] = K[site_xyz,site_xyz] - 1
          if bc_z == 'OBCNeu' and (zp < z or z==0):
            K[site_xyz,site_xyz] = K[site_xyz,site_xyz] - 1
      
          #Coupling terms (don't need to do anything for OBCNeu and OBCDir):
          siteNeigh_x = site(Ltup,xp,y, z )
          siteNeigh_y = site(Ltup,x, yp,z )
          siteNeigh_z = site(Ltup,x, y, zp)
          if (xp > x) or bc_x=='PBC':
            K[siteNeigh_x, site_xyz] = K[siteNeigh_x, site_xyz] - 1.0 
            K[site_xyz, siteNeigh_x] = K[site_xyz, siteNeigh_x] - 1.0
          elif (xp < x) and bc_x=='APBC':
            K[siteNeigh_x, site_xyz] = K[siteNeigh_x, site_xyz] + 1.0 
            K[site_xyz, siteNeigh_x] = K[site_xyz, siteNeigh_x] + 1.0
          
          if (yp > y) or bc_y=='PBC':
            K[siteNeigh_y, site_xyz] = K[siteNeigh_y, site_xyz] - 1.0 
            K[site_xyz, siteNeigh_y] = K[site_xyz, siteNeigh_y] - 1.0
          elif (yp < y) and bc_y=='APBC':
            K[siteNeigh_y, site_xyz] = K[siteNeigh_y, site_xyz] + 1.0 
            K[site_xyz, siteNeigh_y] = K[site_xyz, siteNeigh_y] + 1.0
          
          if (zp > z) or bc_z=='PBC':
            K[siteNeigh_z, site_xyz] = K[siteNeigh_z, site_xyz] - 1.0 
            K[site_xyz, siteNeigh_z] = K[site_xyz, siteNeigh_z] - 1.0
          elif (zp < z) and bc_z=='APBC':
            K[siteNeigh_z, site_xyz] = K[siteNeigh_z, site_xyz] + 1.0 
            K[site_xyz, siteNeigh_z] = K[site_xyz, siteNeigh_z] + 1.0
        #end loop over z
      #end loop over y
    #end loop over x
    
    #use eigh because we know the matrix is symmetric
    Eval,Evec = np.linalg.eigh(K) #, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False, check_finite=True)
    #print Evec
    X = 1./2. * np.matrix(Evec) * np.matrix(np.diag(1. / np.sqrt(Eval))) * np.matrix(Evec.T)
    P = 1./2. * np.matrix(Evec) * np.matrix(np.diag(np.sqrt(Eval))) * np.matrix(Evec.T)
    
    #Make a list of the sites in this region A:
    
    sitesA = np.zeros(LA_x*LA_y*LA_z).tolist() #the indices of the sites in region A
    left_x = (Lx-LA_x)/2
    left_y = (Ly-LA_y)/2
    left_z = (Lz-LA_z)/2
    count = 0
    for x in range(left_x,left_x+LA_x):
      for y in range(left_y,left_y+LA_y):
        for z in range(left_z,left_z+LA_z):
          sitesA[count] = site(Ltup,x,y,z)
          count = count + 1
      
    XA = X[sitesA][:,sitesA]
    PA = P[sitesA][:,sitesA]
  
  else:
    print "*** Boundary condition %sx_%sy_%sz is not supported in getXAPA_edge_K ***" %(bc_x,bc_y,bc_z)
  
  return XA,PA

######################################  getEntropy_Evs  #######################################
###############################################################################################
def getEntropy_Evs(XA,PA):
  #Calculate the matrix CA and its eigenvalues:
  CA_sq = np.matrix(XA)*np.matrix(PA)
  Ev = np.sqrt(np.linalg.eigvals(CA_sq)) #spectrum of eigenvalues of CA_sq
  Ev_new = np.array([e.real for e in Ev if (e.real - 1./2.)>0])

  return Ev_new
  
#  #Calculate the EE for each alpha:
#   S_alpha = np.zeros(len(alpha))
#   for i, n in enumerate(alpha):
#     if n == 1:
#       S_alpha[i] = np.sum( (Ev_new+1./2)*np.log(Ev_new+1./2.) - (Ev_new-1./2.)*np.log(Ev_new-1./2) )
#     else:
#       S_alpha[i] = 1.0/(n-1.0)*np.sum( np.log( (Ev_new+1./2)**n - (Ev_new-1./2.)**n ) )
#   #end alpha loop
#   
#   return S_alpha
#........................................END getEntropy........................................

###############################################################################################
###########################################  site  ############################################
###############################################################################################
def site((Lx,Ly,Lz),x,y,z):
  return x + (y*Lx) + (z*Lx*Ly) #convert (x,y,z) pair to a single site number

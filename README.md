# 2DFEM

P1 FEM method. 
January 2025

Implementation 

1)   FEMmesh is a struct for a P1 mesh that now contains the following fields:  
  
coord     nNodes x 2, coordinates of the nodes  
tr        nTr x 3, counterclockwise-oriented triangles  
eD        neD x 2, Dirichlet edges     
eI        neI x 2, inner edges  
eN        neN x 2, Neumann edges  
eB        neB x 2, boundary edges (not used in this implementation but   
          introduced as a potential replacement for eD and eN)  
normal    neN x 2, non-normalized normal vectors for Neumann edges  
detBk     nTr x 1, determinant of the affine mapping for each triangle  
          (detBk is twice the area of each triangle)      
domain    A domain indicator for each triangle  
domBd     A domain indicator for boundary elements in eB   
b11       \
b12        | entries of the affine mapping matrix   
b21        | (used in mass matrix computation)  
b22       /
c11       \    
c12        | entries of the C_k matrix (used in stress matrix assembly)  
c22       /  
eD2tr     neD x 1, maps the triangle to which the Dirichlet edge belongs     
eN2tr     neN x 1, same as above for Neumann edges       
eB2tr     maps the boundary elements (edges) to their corresponding triangles  
          (intended to replace eD2tr and eN2tr in future updates)   
eI2tr     neI x 2, maps each inner edge to its two associated triangles    
tr2E      nTr x 3, maps edges for each triangle    

Currently, only the following fields are used:   
    coord, tr, eB, eD, eN, normal, detBk, b11, b12, b21, b22, c11, c12, c22  
  
2)   A structure for a P2 mesh, called FEMmeshP2, is very similar to the   
     previous one, with the following differences:   
 
tr        nTr x 6, the first three entries are the vertices; the fourth,  
          fifth, and sixth entries point to the midpoints of the edges opposite  
          to the first, second, and third vertices, respectively.   
eD        neD x 3, Dirichlet edges, with the first two entries as endpoints   
                   and the third as the midpoint  
eI        neI x 3, inner edges with the same structure  
eN        neN x 3, same as above  
eB        neB x 3, same as above   
   
  
2)    Linear FEM functions:   

a)   "FEMStressMatrix"               P1 mass and stress matrices for the Helmholtz   
     "FEMMassMatrix"                 equation.   

b)   "FEMStressMatrixV"              P1 mass and stress matrices for the Helmholtz   
     "FEMMassMatrixV"                equation. VECTORISED VERSION.  
     
c)   "FEMLoadVector"                 P1 load and traction vector for the Helmholtz   
     "FEMTractionVector"             equation.   
     
d)   "FEMLoadVectorV"                P1 load and traction vector for the Helmholtz   
     "FEMTractionVectorV"            equation. VECTORISED VERSION.    

e)   "FEMNonConstantMassMatrix"      P1 mass matrix for a non-constant reaction  
     "FEMNonConstantMassMatrixV"     term  
     
f)   "FEMBoundaryMassMatrix"         P1 boundary mass matrix, suitable for Robin-type  
     "FEMBoundaryMassMatrixV"        boundary conditions (tested but not used here)   

    
3)   Quadratic FEM functions (Prototypes only! Not implemented here)  

a)   "FEMStressMatrixP2"             P2 mass and stress matrices for the Helmholtz 
     "FEMMassMatrixP2"               equation. 

b)   "FEMStressMatrixVP2"            P2 mass and stress matrices for the Helmholtz 
     "FEMMassMatrixVP2"              equation. VECTORISED VERSION.
     
c)   "FEMLoadVectorP2"               P2 load and traction vector for the Helmholtz 
     "FEMTractionVectorP2"           equation. 
     
d)   "FEMLoadVectorVP2"              P2 load and traction vector for the Helmholtz 
     "FEMTractionVectorVP2"          equation. VECTORISED VERSION.
 
4)   For the Stokes problem (P2-P1 Taylor-Hood element)

a)   "FEMmatricesP2P1"    	      Prototypes... 
     "FEMmatricesP2P1V"

4) Meshes

a)   "prepareGridPk"  Constructs a Pk mesh (k = 2,3,4) from a P1 mesh.  
                      Only P2 grids will be used.   
 
b)   "prepareGrid"    Constructs a P1 mesh (FEMmesh) using data obtained   
                      from the MATLAB PDE Toolbox.   

c)   "refineGrid"     Performs RED refinement of a grid (each triangle   
                      is divided into four similar triangles).  
     
5)   Script files (for testing purposes)

 
a)    "scriptTestP1"                 P1 (linear) script with an exact solution for   
                                     testing, including plots, etc. FOR-BASED VERSION.  
                                     
b)    "scriptTestP1V"                P1 (linear) script with an exact solution for   
                                     testing, including plots, etc. VECTORISED VERSION.  
                                     
c)    "scriptTestP2"                 P2 (quadratic), FOR-BASED VERSION.  
                                     
d)    "scriptTestP2V"                P2 (quadratic), VECTORISED VERSION.  

e)    "scriptStokes"                 For the Stokes problem using Taylor-Hood elements   
                                     (P2 for velocity, P1 for pressure).  

Script.  

a)   "PickDomain"     Script used to extract a mesh for a subdomain,   
                      for "anima iocandi". Used for the Pikachu grid   
                      (see below). Nodes on the boundary are always   
                      assumed to be Dirichlet.   

Other utilities  

a)   "LoadMesh"       Script used to switch between different   
                      meshes. It is used in scriptTestP1, scriptTestP1V,   
                      scriptTestP2, and scriptTestP2V.   

Meshes: (All P1)  

a)  "first.mat"       Pure Dirichlet problem  
b)  "second.mat"      Mixed Neumann-Dirichlet problem  
c)  "spanner.mat"     Spanner-like object  
                      [some might say it looks like a rocket ;-) ]  
d)  "mallado6.mat"    Pipe-like structure with a step in the middle.  
                      Suitable for the Stokes problem.  
e)  "SIgrid.mat"      A fun example :-)))  	
g)  "pikachuGrid.mat" A rectangular grid with a hidden Pikachu [ :-O ].  
                      Use PickDomain with indQ = 1 or 2.  
h)  "PuncturedDomain" A domain with a hole. T.eD and T.eN are not defined.  
                      T.eB and T.domBd can be used instead.  
i)  "kinderEgg"       A surprising domain. T.eD and T.eN are not defined.  			
		              T.eB and T.domBd can be used instead.  

Use LoadMesh to load the grid.  

February 2025

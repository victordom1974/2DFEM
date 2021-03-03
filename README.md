# 2DFEM

P1 FEM method. 
January 2021

Implementation 

1)   FEMmesh is a struct for a P1 mesh which has NOW the following fields

coord     nNodes x 2, coordinates of the nodes
tr        nTr x 3, triangles counterclockwise oriented
eD        neD x 2, Dirichlet edges 
eI        neI x 2, Inner edges
eN        neN x 2, Neumann edges
eB        neB x 2  Boundary edges. It's not used in this implementation but 
          is introduced with aim to substitute eD and eN. 
normal    neN x 2, non-normalized normal vectors for Neuamann edges
detBk     nTr x 1, determinant of the affine mapping for each triangle
          (detBk is twice the area of each triangle)    
domain    A domain indicator for each triangle
domBd     A domain indicator for boundary elements in eB 
b11       \
b12        | entries of the matrix of the affine mapping 
b21        | (used in the computation of the mass matrix)
b22       /
c11       \
c12        | entries of the C_k matrix (used in the stress matrix assembly
c22       /
eD2tr     neD x 1 maps the triangle the Dirichlet edge belongs to   
eN2tr     neN x 1 idem for Neuman edges     
eB2tr     maps the boundary elements (edges) to the triangles they belong to
          (is aimed to substitute eD2tr and eN2tr in future updates) 
eI2tr     neI x 2 maps the triangles (2) a inner edges belongs to  
tr2E      nTr x 3 maps the edges for each triangle

Currently, only the following fields are used: 
    coord, tr,eB, eD, eN, normal, detBk, b11, b12, b21, b22, c11, c12, c22

2)   a structure for a P2 mesh, call this FEMmeshP2, is very similar to the 
     former, with these only differences 
 
tr        nTr x 6, the three first entries are the vertices, the fourth
          fith, and sixth pointing to the mid point of each edge opposite to 
          the first, second and third vertex. 
eD        neD x 3, Dirichlet edges, first the end points, the third being
                   the mid point
eI        neI x 3, Inner edges with the same structure
eN        neN x 3, idem
eB        neB x 3, idem 
 


2)    Linear FEM functions: 

a)   "FEMStressMatrix"               P1 Mass and Stress matrices for Helmholtz 
     "FEMMassMatrix"                 equation. 

b)   "FEMStressMatrixV"              P1 Mass and Stress matrices for Helmholtz 
     "FEMMassMatrixV"                equation. VECTORIZED VERSION.
     
c)   "FEMLoadVector"                 P1 load and traction vector for Helmholtz 
     "FEMTractionVector"             equation. 
     
d)   "FEMLoadVectorV"                P1 load and traction vector for Helmholtz 
     "FEMTractionVectorV"            equation. VECTORIZED VERSION.  

e)   "FEMNonConstantMassMatrix"      P1 Mass matrix for a non-constant reaction
     "FEMNonConstantMassMatrixV"     term
     
f)   "FEMBoundaryMassMatrix"         P1 Boundary Mass matrix. Suitable for Robin-type
     "FEMBounaryMassMatrixV"         Boundary conditions (tested, not used here) 

    
3)   Quadratic FEM functions (Just prototipes!!!!)

a)   "FEMStressMatrixP2"             P2 Mass and Stress matrices for Helmholtz 
     "FEMMassMatrixP2"               equation. 

b)   "FEMStressMatrixVP2"            P2 Mass and Stress matrices for Helmholtz 
     "FEMMassMatrixVP2"              equation. VECTORIZED VERSION.
     
c)   "FEMLoadVectorP2"               P2 load and traction vector for Helmholtz 
     "FEMTractionVectorP2"           equation. 
     
d)   "FEMLoadVectorVP2"              P2 load and traction vector for Helmholtz 
     "FEMTractionVectorVP2"          equation. VECTORIZED VERSION.
 
4)   For Stokes problem (P2-P1 Taylor hood element)

a)   "FEMmatricesP2P1"    	      Prototipes... 
     "FEMmatricesP2P1V"

4) Meshes

a)   "prepareGridPk"  construct, from a P1 mesh, a Pk mesh, with k =2,3,4. 
                      We will use only the P2 grids   
 
b)   "prepareGrid"    construct a P1 mesh (FEMmesh) from data obtained 
                      from Matlab PDETOOL toolbox 

c)   "refineGrid"     preforms a RED refinement of a grid (four new similar
                      triangles for each one) 
     
5)   Scritp files (for testing purposes)

 
a)    "scriptTestP1"                 P1 (linear) script with exact solution to test with, 
                                     plots, and so on. FOR-BASED VERSION. 
                                     
b)    "scriptTestP1V"                P1 (linear) script with exact solution to test with, 
                                     plots, and so on. VECTORIZED VERSION.
                                     
c)    "scriptTestP2V"                P2 (quadratic), FOR-BASED VERSION. 
                                     
d)    "scriptTestP2V"                P2 (quadratic), VECTORIZED VERSION.  

e)    "scriptStokes"                 For stokes problem, Taylor-Hoold elements, 
                                     P2 for velocity, P1 for pressure 
                                     

Script.  


a)   "PickDomain"     script file used to take a mesh for a subdomain, for 
                      "anima iocandi". Used for pikachu grid (see below). 
                      Nodes on the boundary are assumed to be Dirichlet 
                      in all the cases. 

Some uncategorized stuff

a)   "PickDomain"     script file used to take a mesh for a subdomain, for 
                      "anima iocandi". Used for pikachu grid (see below). 
                      Nodes on the boundary are assumed to be Dirichlet 
                      in all the cases.                       
                      
a)   "LoadMesh"       script file used to switch among the different
                      meshes. It is used in scriptTestP1, scriptTestP1V,
                      scriptTestP2 and scriptTestP2V files. 

Meshes: (P1 all of them)

a)  "first.mat"       a pure Dirichlet problem  
b)  "second.mat"      mixed Neumann-Dirichlet problem 
c)  "spanner.mat"     a spanner-like object
                      [some might say this is more like a rocket ;-) ]
d)  "mallado6.mat"    a sort of pipe with a step in the middle. Suitable for Stokes
                      problem, for example. 
e)  "SIgrid.mat"      a funny example :-)))
f)  "Elephant.mat"    well, what can it be? an elephant, of course :-)))))
g)  "pikachuGrid.mat" A rectangular grid with a hidden Pikachu [ :-O ]. 
                      Use PickDomain with indQ =1 or 2.
h)  "PuncturedDomain" just a domain with a hole. T.eD and T.eN are not defined.
                      T.eB and T.domBd can be used to set as desired. 
i)  "kinderEgg"       Just a surprising domain it is.  T.eD and T.eN are not defined.
		              T.eB and T.domBd can be used instead.


Use Load ***** to load the grid


March 2021

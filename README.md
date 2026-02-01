# 2DFEM

P1 finite element method on 2D triangular meshes.

Author: Victor Dominguez Baguena  
Date: January 2026

--------------------------------------------------------------------------

1. Mesh data structure (P1)

The finite element mesh is stored in a struct called FEMmesh, describing a
P1 (linear) triangular mesh with the following fields.

Geometry and connectivity:
- coord   : nNodes x 2 array, nodal coordinates.
- tr      : nTr x 3 array, triangle connectivity (counterclockwise oriented).

Boundary and interior edges:
- eD      : neD x 2 array, Dirichlet edges.
- eN      : neN x 2 array, Neumann edges.
- eI      : neI x 2 array, interior edges.
- eB      : neB x 2 array, boundary edges (general boundary set).
            This field is intended to replace eD and eN in future updates.
- normal  : neN x 2 array, outward normals on Neumann edges, scaled by the
            edge length (non-unit normals).

Element geometry and mappings:
- detBk   : nTr x 1 array, determinant of the affine mapping for each triangle.
            detBk equals twice the area of the triangle.
- domain  : nTr x 1 array, domain indicator for each triangle.
- domB    : neB x 1 array, domain indicator for boundary edges in eB.

Affine-map related coefficients:
- b11, b12, b21, b22 : entries of the affine mapping matrix B_K.
- c11, c12, c22      : entries of the matrix C_K used in the stiffness
                       (stress) matrix assembly.

Edge-to-triangle and triangle-to-edge maps:
- eD2tr   : neD x 1 array, triangle adjacent to each Dirichlet edge.
- eN2tr   : neN x 1 array, triangle adjacent to each Neumann edge.
- eB2tr   : boundary-edge to triangle map (intended to replace eD2tr/eN2tr).
- eI2tr   : neI x 2 array, triangles sharing each interior edge.
- tr2E    : nTr x 3 array, edge indices associated with each triangle.

Currently used fields in the P1 Helmholtz implementation:
coord, tr, eB, eD, eN, normal, detBk,
b11, b12, b21, b22, c11, c12, c22.

--------------------------------------------------------------------------

2. P1 FEM routines for the Helmholtz equation

Stiffness and mass matrices:
- FEMStressMatrix
- FEMMassMatrix

Vectorized versions:
- FEMStressMatrixV
- FEMMassMatrixV

Right-hand side assembly:
- FEMLoadVector
- FEMTractionVector

Vectorized versions:
- FEMLoadVectorV
- FEMTractionVectorV

Non-constant reaction term:
- FEMNonConstantMassMatrix
- FEMNonConstantMassMatrixV

Boundary mass matrix (Robin-type boundary terms):
- FEMMassBoundaryMatrix
- FEMMassBoundaryMatrixV

--------------------------------------------------------------------------

3. Mesh utilities

- prepareGrid   : construct a P1 mesh (FEMmesh) from MATLAB PDETOOL data.
- refineGrid    : uniform RED refinement (each triangle is subdivided into
                  four similar triangles).

--------------------------------------------------------------------------

4. Test scripts

- scriptTestP1   : P1 test with exact solution, plots and error computation.
                   Assembly performed with for-loops.
- scriptTestP1V  : Same test using vectorized assembly routines.

--------------------------------------------------------------------------

5. Additional scripts

- PickDomain : extract a subdomain from a mesh. Boundary nodes are assumed
               to be Dirichlet.
- loadMesh   : script used to switch among different predefined meshes.
               Used in scriptTestP1 and scriptTestP1V.

--------------------------------------------------------------------------

6. Available meshes (P1)

- first.mat        : pure Dirichlet problem.
- second.mat       : mixed Dirichlet–Neumann problem.
- spanner.mat      : spanner-like geometry.
- mallado6.mat     : pipe-like domain, suitable for Stokes-type tests.
- SIgrid.mat       : example geometry.
- Elephant.mat     : elephant-shaped domain.
- pikachuGrid.mat  : rectangular grid with a hidden Pikachu.
                     Use PickDomain with indQ = 1 or 2.
- PuncturedDomain  : domain with a hole. eD and eN may be missing; eB and domB
                     can be used instead.
- kinderEgg        : domain with similar setup; eD and eN may be missing;
                     eB and domB can be used instead.


% MESH  Triangular, tetrahedral, polygon, and surface meshes with marking strategies
%
%   Classes
%   Mark               -  Base class of mesh marking strategies for adaptive refinement
%   MarkBulk           -  Bulk marking strategy with parameter theta
%   MarkFixedFraction  -  Fixed fraction marking strategy scaled by the polynomial degree
%   MarkHP             -  NVB based hp marking strategy bisecting one edge
%   MarkHP3Edges       -  NVB based hp marking strategy bisecting three edges
%   MarkMaximum        -  Maximum marking strategy with parameter theta
%   MarkMean           -  Mean marking strategy of the indicators
%   MarkP              -  Pure p-refinement marking strategy
%   MarkUniform        -  Uniform marking strategy
%   Mesh               -  Triangular mesh with adaptive refinement and coarsening
%   Mesh3d             -  Tetrahedral volume mesh with recursive newest vertex bisection
%   MeshPoly           -  Polygon mesh with arbitrary polygonal elements
%   MeshPolyH          -  Polygon mesh with hanging-node-free H-refinement
%   MeshPolyQ          -  Quadrilateral mesh stored as polygons with a quad mesh
%   MeshPolyT          -  Triangular mesh stored as polygons with a triangulation
%   MeshSurface        -  Triangular surface mesh embedded in 3d
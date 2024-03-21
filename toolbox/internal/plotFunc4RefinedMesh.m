function plotFunc4RefinedMesh(f,mesh,sub)
%plotFunc4RefinedMesh - Plot element function on subtriangulation
%   plotFunc4RefinedMesh(F,MESH) plots the function F over the 2d
%   triangular mesh MESH with patch, one patch object per triangle of
%   the subtriangulation with one face per element, with the values of F
%   at the vertices as height and color. F(X) receives an nrElems-by-2
%   array X whose row i is a point of element i of MESH and returns the
%   nrElems-by-1 values at these points; F must be scalar valued.
%
%   plotFunc4RefinedMesh(F,MESH,SUB) refines the subtriangulation SUB
%   times uniformly (default 0, i.e. one triangle per element).
%
%   See also plotRefFunc4RefinedMesh, patch

% Copyright (c) 2026 by the fem4m authors

if nargin < 3 || isempty(sub)
    sub = 0;
end
meshRef = Mesh([0,0;1,0;0,1],[2 3 1],...
    [1 2; 2 3; 3 1],[]);
if sub > 0
    meshRef = meshRef.refine('uniform',sub);
end
z = zeros(mesh.nrElems,meshRef.nrNodes);
c1 = mesh.c4n(mesh.n4e(:,1),:);
c2 = mesh.c4n(mesh.n4e(:,2),:);
c3 = mesh.c4n(mesh.n4e(:,3),:);
for j = 1: meshRef.nrNodes
    z(:,j) = f(c1 + meshRef.c4n(j,1)*(c2-c1) + meshRef.c4n(j,2)*(c3-c1));
end

coordXref = reshape(meshRef.c4n(meshRef.n4e,1),[],3);
coordYref = reshape(meshRef.c4n(meshRef.n4e,2),[],3);
coordX = reshape(mesh.c4n(mesh.n4e,1),[],3);
coordY = reshape(mesh.c4n(mesh.n4e,2),[],3);

for j = 1 : meshRef.nrElems
    X(1,:) = coordX(:,1) + (coordX(:,2)-coordX(:,1))*coordXref(j,1)+(coordX(:,3)-coordX(:,1))*coordYref(j,1);
    X(2,:) = coordX(:,1) + (coordX(:,2)-coordX(:,1))*coordXref(j,2)+(coordX(:,3)-coordX(:,1))*coordYref(j,2);
    X(3,:) = coordX(:,1) + (coordX(:,2)-coordX(:,1))*coordXref(j,3)+(coordX(:,3)-coordX(:,1))*coordYref(j,3);
    Y(1,:) = coordY(:,1) + (coordY(:,2)-coordY(:,1))*coordXref(j,1)+(coordY(:,3)-coordY(:,1))*coordYref(j,1);
    Y(2,:) = coordY(:,1) + (coordY(:,2)-coordY(:,1))*coordXref(j,2)+(coordY(:,3)-coordY(:,1))*coordYref(j,2);
    Y(3,:) = coordY(:,1) + (coordY(:,2)-coordY(:,1))*coordXref(j,3)+(coordY(:,3)-coordY(:,1))*coordYref(j,3);
    curZ = z(:,meshRef.n4e(j,:));
    patch(X,Y,curZ',curZ')
end
end
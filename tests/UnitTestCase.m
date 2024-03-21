classdef (Abstract) UnitTestCase < PatchTestCase
    %UnitTestCase - Common base class of the fem4m unit tests
    %   UnitTestCase extends PatchTestCase, which adds the fem4m toolbox to
    %   the path, by the reference elements, meshes and polynomials shared
    %   by the unit and convergence tests of the toolbox building blocks.
    %
    %   UnitTestCase functions:
    %       referenceTriangle    - Mesh of the reference triangle
    %       referenceTetrahedron - Mesh of the reference tetrahedron
    %       cubeMesh             - Unit cube mesh refined by bisection
    %       polynomial           - Polynomial of total degree k in 2d
    %
    %   See also PatchTestCase

    % Copyright (c) 2026 by the fem4m authors

    methods (Static)

        function mesh = referenceTriangle()
            %referenceTriangle - Mesh of the reference triangle
            %   MESH = UnitTestCase.referenceTriangle() is the one-element
            %   mesh of conv{(0,0),(1,0),(0,1)} with Dirichlet boundary.
            mesh = Mesh([0 0;1 0;0 1],[1 2 3],[1 2;2 3;3 1],[]);
        end

        function mesh = referenceTetrahedron()
            %referenceTetrahedron - Mesh of the reference tetrahedron
            %   MESH = UnitTestCase.referenceTetrahedron() is the
            %   one-element mesh of conv{0,e1,e2,e3} with Dirichlet
            %   boundary.
            mesh = Mesh3d([0 0 0;1 0 0;0 1 0;0 0 1],[1 2 3 4], ...
                [1 3 2;1 2 4;1 4 3;2 3 4],[],[],'t4e',0);
        end

        function mesh = cubeMesh(bdry)
            %cubeMesh - Unit cube mesh refined by bisection
            %   MESH = UnitTestCase.cubeMesh() bisects all elements of the
            %   unit cube mesh once and then the first element again, so
            %   that the mesh is not uniform.
            %
            %   MESH = UnitTestCase.cubeMesh(BDRY) also sets the boundary
            %   conditions BDRY of the six faces, e.g. 'DNDNDN'.
            opts = struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1);
            if nargin > 0
                opts.bdry = bdry;
            end
            mesh = Mesh3d(Geometry(opts));
            mesh = mesh.refine(true(mesh.nrElems,1));
            mesh = mesh.refine([true;false(mesh.nrElems-1,1)]);
        end

        function f = polynomial(k)
            %polynomial - Polynomial of total degree k in two variables
            %   F = UnitTestCase.polynomial(K) is the function handle of
            %   x^k - 2*x*y^(k-1) + y + 1 for points given as rows.
            f = @(z)z(:,1).^k - 2*z(:,1).*z(:,2).^(k-1) + z(:,2) + 1;
        end

    end

end
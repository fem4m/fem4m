classdef (Abstract) PatchTestCase < matlab.unittest.TestCase
    %PatchTestCase - Common base class of the fem4m patch tests
    %   PatchTestCase adds the fem4m toolbox to the path and provides the
    %   small, locally refined test meshes used by the patch tests. A patch
    %   test solves a problem whose polynomial exact solution lies in the
    %   discrete space and checks that it is reproduced up to round-off.
    %
    %   PatchTestCase properties:
    %       tol         - Absolute tolerance for errors of direct solvers
    %
    %   PatchTestCase functions:
    %       patchMesh   - Uniformly and locally refined 2d test mesh
    %
    %   See also matlab.unittest.TestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (Constant)
        tol = 1e-9
    end

    methods (TestClassSetup)

        function addToolboxPath(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,'toolbox'),'IncludingSubfolders',true));
        end

    end

    methods (Static)

        function mesh = patchMesh(opts,nrUniform)
            %patchMesh - Uniformly and locally refined 2d test mesh
            %   MESH = PatchTestCase.patchMesh(OPTS) creates the mesh of the
            %   geometry Geometry(OPTS), refines it uniformly once and then
            %   refines the first element locally by red-green-blue
            %   refinement, so that the mesh is not uniform.
            %
            %   MESH = PatchTestCase.patchMesh(OPTS,NRUNIFORM) refines
            %   NRUNIFORM times uniformly before the local refinement.
            if nargin < 2
                nrUniform = 1;
            end
            mesh = Mesh(Geometry(opts));
            if nrUniform > 0
                mesh = mesh.refine('uniform',nrUniform);
            end
            marked4e = zeros(mesh.nrElems,1);
            marked4e(1) = 1;
            mesh = mesh.refine(marked4e);
        end

    end

end
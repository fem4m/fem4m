classdef exLaplaceBeltramiPkAdaptivePatchTest < PatchTestCase
    %exLaplaceBeltramiPkAdaptivePatchTest - Patch test of the example exLaplaceBeltramiPkAdaptive
    %   Conforming Pk finite elements for the Laplace-Beltrami operator on
    %   the flat surface mesh of problemSurfaceSquareExact reproduce its
    %   biquadratic solution for k >= 4; the residual estimator and the
    %   geometric oscillation vanish. On curved surfaces the geometric
    %   error prevents an exact reproduction.
    %
    %   See also exLaplaceBeltramiPkAdaptive, MeshSurface,
    %   SurfaceTransformation, problemSurfaceSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degree = struct('P4',4,'P5',5)
    end

    properties
        mesh
        p
    end

    methods (TestClassSetup)

        function createMesh(testCase)
            testCase.p = problemSurfaceSquareExact();
            p = testCase.p;
            testCase.mesh = MeshSurface(p.geom,@p.surfaceFunction,@p.gradSurfaceFunction,1);
            testCase.mesh = testCase.mesh.refine('uniform');
            marked4e = zeros(testCase.mesh.nrElems,1);
            marked4e(1) = 1;
            testCase.mesh = testCase.mesh.refine(marked4e);
        end

    end

    methods (Test)

        function testL2Error(testCase,degree)
            p = testCase.p;
            [x,fe] = exLaplaceBeltramiPkAdaptivePatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,2*(degree+1),'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testH1Error(testCase,degree)
            p = testCase.p;
            [x,fe] = exLaplaceBeltramiPkAdaptivePatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,2*(degree+1),'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testEstimator(testCase,degree)
            p = testCase.p;
            [x,fe,q] = exLaplaceBeltramiPkAdaptivePatchTest.solve(testCase.mesh,p,degree);

            eta = EstimatorForm(['h^2*(f(x)+div(grad(u)))^2dx', ...
                '+ h^1*([grad(u)*n])^2dsi']);
            eta4e = eta.assemble(x,fe,'f',@(pts,pts_ref,parts)(q.f(pts,pts_ref,parts)));

            testCase.verifyLessThan(norm(eta4e),testCase.tol);
        end

        function testGeometricOscillation(testCase,degree)
            p = testCase.p;
            [x,fe,q] = exLaplaceBeltramiPkAdaptivePatchTest.solve(testCase.mesh,p,degree);

            mu = EstimatorForm('(B(x)*grad(u))^2dx');
            mu4e = mu.assemble(x,fe,'B',@(pts,pts_ref,parts)(q.B(pts,pts_ref,parts)));

            testCase.verifyLessThan(norm(mu4e),testCase.tol);
        end

    end

    methods (Static)

        function [x,fe,q] = solve(mesh,p,degree)
            fe = Pk(degree,mesh);
            q = SurfaceTransformation(mesh,p);
            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);
            b = LinearForm('f(x) v dx').assemble(fe,'f',@(pts,pts_ref,parts)(q.f(pts,pts_ref,parts)));
            x = zeros(size(A,1),1);
            x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
            b = b - A*x;
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
        end

    end

end
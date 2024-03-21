classdef exPoisson3dPkAdaptivePatchTest < PatchTestCase
    %exPoisson3dPkAdaptivePatchTest - Patch test of the example exPoisson3dPkAdaptive
    %   Conforming Pk finite elements on tetrahedra with inhomogeneous
    %   Dirichlet and Neumann data reproduce every monomial
    %   u = x^i y^j z^l with i+j+l <= k, and the residual-based error
    %   estimator vanishes.
    %
    %   See also exPoisson3dPkAdaptive, problemScalarEllipticCubeMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        monom = struct( ...
            'P1_x1',struct('degree',1,'i',1,'j',0,'l',0), ...
            'P1_z1',struct('degree',1,'i',0,'j',0,'l',1), ...
            'P2_x1z1',struct('degree',2,'i',1,'j',0,'l',1), ...
            'P2_y2',struct('degree',2,'i',0,'j',2,'l',0), ...
            'P3_x1y1z1',struct('degree',3,'i',1,'j',1,'l',1), ...
            'P3_z3',struct('degree',3,'i',0,'j',0,'l',3))
    end

    properties
        mesh
    end

    methods (TestClassSetup)

        function createMesh(testCase)
            geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5], ...
                'length',1,'bdry','DNDNDN'));
            testCase.mesh = Mesh3d(geom);
            testCase.mesh = testCase.mesh.refine(true(testCase.mesh.nrElems,1));
            testCase.mesh = testCase.mesh.refine(true(testCase.mesh.nrElems,1));
            testCase.mesh = testCase.mesh.refine([true;false(testCase.mesh.nrElems-1,1)]);
        end

    end

    methods (Test)

        function testH1Error(testCase,monom)
            p = problemScalarEllipticCubeMonomExact(monom.i,monom.j,monom.l);
            [x,fe] = exPoisson3dPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testL2Error(testCase,monom)
            p = problemScalarEllipticCubeMonomExact(monom.i,monom.j,monom.l);
            [x,fe] = exPoisson3dPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testEstimator(testCase,monom)
            p = problemScalarEllipticCubeMonomExact(monom.i,monom.j,monom.l);
            [x,fe] = exPoisson3dPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            eta = EstimatorForm(['h^2*(f(x)+div(grad(u)))^2dx', ...
                '+ h^1*([grad(u)*n])^2dsi', ...
                '+ h^1*(g(x)-grad(u)*n)^2dsn']);
            eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n));

            testCase.verifyLessThan(norm(eta4e),testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(mesh,p,degree)
            fe = Pk3d(degree,mesh);
            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);
            b = LinearForm('f(x) v dx + g(x) v dsn').assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
            x = zeros(size(A,1),1);
            x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
            b = b - A*x;
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
        end

    end

end
classdef exPoissonHPAdaptivePatchTest < PatchTestCase
    %exPoissonHPAdaptivePatchTest - Patch test of the example exPoissonHPAdaptive
    %   hp finite elements HPk with a varying polynomial degree per element
    %   and inhomogeneous Dirichlet and Neumann data reproduce every
    %   monomial u = x^i y^j whose degree does not exceed the minimal
    %   element degree, and the residual-based error estimator vanishes.
    %
    %   See also exPoissonHPAdaptive, HPk,
    %   problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (ClassSetupParameter)
        domain = struct('Square',struct('name','Square','bdry','DNDN'), ...
            'Lshape',struct('name','Lshape','bdry','DN'))
    end

    properties (TestParameter)
        monom = struct( ...
            'P1P2_x1',struct('minDegree',1,'i',1,'j',0), ...
            'P2P3_x1y1',struct('minDegree',2,'i',1,'j',1), ...
            'P2P3_y2',struct('minDegree',2,'i',0,'j',2), ...
            'P3P4_x2y1',struct('minDegree',3,'i',2,'j',1))
    end

    properties
        mesh
    end

    methods (TestClassSetup)

        function createMesh(testCase,domain)
            testCase.mesh = PatchTestCase.patchMesh(domain);
        end

    end

    methods (Test)

        function testH1Error(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,fe] = exPoissonHPAdaptivePatchTest.solve(testCase.mesh,p,monom.minDegree);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,12,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testEstimator(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,fe] = exPoissonHPAdaptivePatchTest.solve(testCase.mesh,p,monom.minDegree);

            eta = EstimatorForm(['h^2*p^-2*(f(x)+div(grad(u)))^2dx', ...
                '+ h^1*p^-1*([grad(u)*n])^2dsi', ...
                '+ h^1*p^-1*(g(x)-grad(u)*n)^2dsn']);
            eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n));

            testCase.verifyLessThan(norm(eta4e),testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(mesh,p,minDegree)
            % alternate the element degrees between minDegree and minDegree+1
            degree4e = minDegree + mod((1:mesh.nrElems)',2);
            fe = HPk(degree4e,mesh);
            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);
            b = LinearForm('f(x) v dx + g(x) v dsn').assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
            x = fe.projectDbBoundary(@p.u4Db);
            b = b - A*x;
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
        end

    end

end
classdef exPoissonPkAdaptivePatchTest < PatchTestCase
    %exPoissonPkAdaptivePatchTest - Patch test of the example exPoissonPkAdaptive
    %   Conforming Pk finite elements with inhomogeneous Dirichlet and
    %   Neumann data reproduce every monomial u = x^i y^j with i+j <= k,
    %   and the residual-based error estimator vanishes.
    %
    %   See also exPoissonPkAdaptive, problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (ClassSetupParameter)
        domain = struct('Square',struct('name','Square','bdry','DNDN'), ...
            'Lshape',struct('name','Lshape','bdry','DN'))
    end

    properties (TestParameter)
        monom = struct( ...
            'P1_x1',struct('degree',1,'i',1,'j',0), ...
            'P1_y1',struct('degree',1,'i',0,'j',1), ...
            'P2_x1y1',struct('degree',2,'i',1,'j',1), ...
            'P2_y2',struct('degree',2,'i',0,'j',2), ...
            'P3_x2y1',struct('degree',3,'i',2,'j',1), ...
            'P3_x3',struct('degree',3,'i',3,'j',0), ...
            'P4_x2y2',struct('degree',4,'i',2,'j',2), ...
            'P4_x1y3',struct('degree',4,'i',1,'j',3))
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
            [x,fe] = exPoissonPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testL2Error(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,fe] = exPoissonPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testEstimator(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,fe] = exPoissonPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            eta = EstimatorForm(['h^2*p^-2*(f(x)+div(grad(u)))^2dx', ...
                '+ h^1*p^-1*([grad(u)*n])^2dsi', ...
                '+ h^1*p^-1*(g(x)-grad(u)*n)^2dsn']);
            eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n));

            testCase.verifyLessThan(norm(eta4e),testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(mesh,p,degree)
            fe = Pk(degree,mesh);
            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);
            b = LinearForm('f(x) v dx + g(x) v dsn').assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
            x = zeros(size(A,1),1);
            x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
            b = b - A*x;
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
        end

    end

end
classdef exHelmholtzPkAdaptivePatchTest < PatchTestCase
    %exHelmholtzPkAdaptivePatchTest - Patch test of the example exHelmholtzPkAdaptive
    %   Conforming Pk finite elements for the Helmholtz equation
    %   -Delta u - omega^2 u = f with inhomogeneous Dirichlet and impedance
    %   (Robin) data reproduce every monomial u = x^i y^j with i+j <= k,
    %   and the residual-based error estimator vanishes.
    %
    %   See also exHelmholtzPkAdaptive, problemHelmholtzSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        monom = struct( ...
            'P1_x1',struct('degree',1,'i',1,'j',0), ...
            'P2_x1y1',struct('degree',2,'i',1,'j',1), ...
            'P2_y2',struct('degree',2,'i',0,'j',2), ...
            'P3_x2y1',struct('degree',3,'i',2,'j',1), ...
            'P3_x3',struct('degree',3,'i',3,'j',0))
        wavenumber = struct('k1',1,'k3',3)
    end

    properties
        mesh
    end

    methods (TestClassSetup)

        function createMesh(testCase)
            testCase.mesh = PatchTestCase.patchMesh(struct('name','Square','bdry','DNDN'));
        end

    end

    methods (Test)

        function testL2Error(testCase,monom,wavenumber)
            p = problemHelmholtzSquareMonomExact(wavenumber,monom.i,monom.j);
            [x,fe] = exHelmholtzPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,2*(monom.degree+3),'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testH1Error(testCase,monom,wavenumber)
            p = problemHelmholtzSquareMonomExact(wavenumber,monom.i,monom.j);
            [x,fe] = exHelmholtzPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,2*(monom.degree+3),'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testEstimator(testCase,monom,wavenumber)
            p = problemHelmholtzSquareMonomExact(wavenumber,monom.i,monom.j);
            [x,fe] = exHelmholtzPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            omega = p.omega;
            eta = EstimatorForm(['h^2*(f(x)+div(grad(u))+omega2(x)*u)^2dx', ...
                '+ h^1*([grad(u)*n])^2dsi', ...
                '+ h^1*(g(x)-grad(u)*n-iomega(x)*u)^2dsn']);
            eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g4Rb(x,n),'omega2',omega^2,'iomega',1i*omega);

            testCase.verifyLessThan(norm(eta4e),testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(mesh,p,degree)
            omega = p.omega;
            fe = Pk(degree,mesh);
            A = BilinearForm('grad(u)*grad(v)dx - omega2(x)*u*v dx + iomega(x)*u*v dsn');
            A = A.assemble(fe,'omega2',omega^2,'iomega',1i*omega);
            b = LinearForm('f(x) v dx + g(x) v dsn').assemble(fe,'f',@p.f,'g',@(x,n)p.g4Rb(x,n));
            x = zeros(size(A,1),1);
            x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
            b = b - A*x;
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
        end

    end

end
classdef exPoissonDGAdaptivePatchTest < PatchTestCase
    %exPoissonDGAdaptivePatchTest - Patch test of the example exPoissonDGAdaptive
    %   The symmetric interior penalty dG method DGPk with homogeneous
    %   Dirichlet data on the re-entrant edges of the L-shape and
    %   inhomogeneous Neumann data elsewhere reproduces every monomial
    %   u = x^i y^j with i,j >= 1 (which vanishes on the re-entrant edges)
    %   and i+j <= k, and the error estimator vanishes.
    %
    %   See also exPoissonDGAdaptive, DGPk,
    %   problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        monom = struct( ...
            'DGP2_x1y1',struct('degree',2,'i',1,'j',1), ...
            'DGP3_x1y1',struct('degree',3,'i',1,'j',1), ...
            'DGP3_x2y1',struct('degree',3,'i',2,'j',1), ...
            'DGP3_x1y2',struct('degree',3,'i',1,'j',2))
    end

    properties
        mesh
    end

    methods (TestClassSetup)

        function createMesh(testCase)
            testCase.mesh = PatchTestCase.patchMesh(struct('name','Lshape','bdry','DN'));
        end

    end

    methods (Test)

        function testH1Error(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,fe] = exPoissonDGAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('(exactGrad(x)-grad(u))^2dx').assemble(x,fe,2*fe.degree,'exactGrad',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testL2Error(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,fe] = exPoissonDGAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('(exact(x)-u)^2dx').assemble(x,fe,2*fe.degree+2,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testEstimator(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,fe,penalty] = exPoissonDGAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            eta = EstimatorForm(['h^2*(f(x)+div(grad(u)))^2dx', ...
                '+ h^1*([grad(u)*n])^2dsi + h^-1*([alpha(x)*u])^2dsi', ...
                '+ h^1*(g(x)-grad(u)*n)^2dsn + h^-1*([alpha(x)*u])^2dsd']);
            eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n),'alpha',sqrt(penalty));

            testCase.verifyLessThan(norm(eta4e),testCase.tol);
        end

    end

    methods (Static)

        function [x,fe,penalty] = solve(mesh,p,degree)
            penalty = 2*(degree+1)*degree;
            fe = DGPk(degree,mesh);
            a = BilinearForm(['grad(u)*grad(v)dx', ...
                '- [u]*n.{grad(v)}dsi - {grad(u)}.[v]*ndsi + h^-1*alpha(x)*[u][v] dsi', ...
                '- [u]*n.{grad(v)}dsd - {grad(u)}.[v]*ndsd + h^-1*alpha(x)*[u][v] dsd']);
            A = a.assemble(fe,'alpha',penalty);
            b = LinearForm('f(x)vdx + g(x)vdsn').assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
            x = A\b;
        end

    end

end
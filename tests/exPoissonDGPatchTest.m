classdef exPoissonDGPatchTest < PatchTestCase
    %exPoissonDGPatchTest - Patch test of the example exPoissonDG
    %   The symmetric interior penalty dG method DGPk with homogeneous
    %   Dirichlet data imposed by Nitsche's method on all boundary edges
    %   reproduces the biquadratic solution of
    %   problemScalarEllipticSquareExact for k >= 4.
    %
    %   See also exPoissonDG, DGPk, problemScalarEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degree = struct('DGP4',4,'DGP5',5)
    end

    properties
        mesh
    end

    methods (TestClassSetup)

        function createMesh(testCase)
            testCase.mesh = PatchTestCase.patchMesh(struct('name','Square'));
        end

    end

    methods (Test)

        function testL2Error(testCase,degree)
            p = problemScalarEllipticSquareExact();
            [x,fe] = exPoissonDGPatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('(exact(x)-u)^2dx').assemble(x,fe,2*fe.degree+2,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testDGNormError(testCase,degree)
            p = problemScalarEllipticSquareExact();
            [x,fe] = exPoissonDGPatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('(exactGrad(x)-grad(u))^2dx + h^-1*([u])^2ds').assemble(x,fe, ...
                2*fe.degree,'exactGrad',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(mesh,p,degree)
            penalty = 2*(degree+1)*degree;
            fe = DGPk(degree,mesh);
            a = BilinearForm(['grad(u)*grad(v)dx', ...
                '- [u]*n.{grad(v)}ds - {grad(u)}.[v]*nds', ...
                '+ h^-1*alpha(x)*[u][v] ds']);
            A = a.assemble(fe,'alpha',penalty);
            b = LinearForm('f(x)vdx').assemble(fe,'f',@p.f);
            x = A\b;
        end

    end

end
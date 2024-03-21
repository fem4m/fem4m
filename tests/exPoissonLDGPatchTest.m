classdef exPoissonLDGPatchTest < PatchTestCase
    %exPoissonLDGPatchTest - Patch test of the example exPoissonLDG
    %   The local dG method with the discrete gradient G mapping DGDk into
    %   DGDkVector(k-1) and homogeneous Dirichlet data reproduces the
    %   biquadratic solution of problemScalarEllipticSquareExact and its
    %   gradient for k >= 4.
    %
    %   See also exPoissonLDG, problemScalarEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degree = struct('DGD4',4,'DGD5',5)
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
            [x,fe] = exPoissonLDGPatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('(exact(x)-u)^2dx').assemble(x,fe,2*fe.degree+2,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testDGNormError(testCase,degree)
            p = problemScalarEllipticSquareExact();
            [x,fe] = exPoissonLDGPatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('(exactGrad(x)-grad(u))^2 dx + h^-1*([u])^2 ds').assemble(x,fe, ...
                2*fe.degree,'exactGrad',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testDiscreteGradientError(testCase,degree)
            p = problemScalarEllipticSquareExact();
            [x,~,fe2,G] = exPoissonLDGPatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('(exactGrad(x)-u)^2dx').assemble(G*x,fe2,2*fe2.degree, ...
                'exactGrad',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe,fe2,G] = solve(mesh,p,degree)
            eta = 1;
            fe2 = DGDkVector(degree-1,mesh);
            fe = DGDk(degree,mesh);
            A = BilinearForm('sigma*tau dx',{'sigma','tau'}).assemble(fe2);
            B = BilinearForm('grad(u)*tau dx - [u]*n{tau} ds ',{'u','tau'}).assemble(fe,fe2);
            G = A\B;
            S = BilinearForm('h^-1*eta(x)*[u][v] ds').assemble(fe,'eta',eta);
            b = LinearForm('f(x)vdx').assemble(fe,'f',@p.f);
            x = (G'*A*G + S)\b;
        end

    end

end
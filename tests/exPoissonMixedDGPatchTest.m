classdef exPoissonMixedDGPatchTest < PatchTestCase
    %exPoissonMixedDGPatchTest - Patch test of the example exPoissonMixedDG
    %   The mixed dG method DGDkVector x DGDk with the numerical flux
    %   parameters C11 = C22 = 1, C12 = 0 and homogeneous Dirichlet data
    %   reproduces the biquadratic solution of
    %   problemScalarEllipticSquareExact and sigma = -grad(u) for k >= 4.
    %
    %   See also exPoissonMixedDG, problemScalarEllipticSquareExact

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
            [x,fe,dg] = exPoissonMixedDGPatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('(exact(x)-u)^2dx').assemble(x(fe.nrDoF+1:end),dg, ...
                2*fe.degree+2,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testSigmaError(testCase,degree)
            p = problemScalarEllipticSquareExact();
            [x,fe] = exPoissonMixedDGPatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('(exactGrad(x)+sigma)^2dx + (C(x)[sigma])^2dsi','sigma').assemble( ...
                x(1:fe.nrDoF),fe,2*fe.degree+2,'exactGrad',@p.exactGradU,'C',1);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testJumpError(testCase,degree)
            p = problemScalarEllipticSquareExact();
            [x,fe,dg] = exPoissonMixedDGPatchTest.solve(testCase.mesh,p,degree);

            err = ErrorForm('(C(x)[u])^2ds').assemble(x(fe.nrDoF+1:end),dg,2*fe.degree,'C',1);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe,dg] = solve(mesh,p,degree)
            C11 = 1;
            C12 = [0,0];
            C22 = 1;
            fe = DGDkVector(degree,mesh);
            dg = DGDk(degree,mesh);
            A = BilinearForm('sigma*tau dx + C(x)*[sigma]*n [tau]*n dsi',{'sigma','tau'}).assemble(fe,'C',C22);
            B = BilinearForm('div(sigma) v dx - [sigma]*n{v} dsi',{'sigma','v'}).assemble(fe,dg);
            B = B + BilinearForm('C(x)[u]*n [tau]*n dsi',{'u','tau'}).assemble(dg,fe,'C',C12)';
            S = BilinearForm('C(x)*[u][v] ds').assemble(dg,'C',C11);
            b = LinearForm('f(x)vdx').assemble(dg,'f',@p.f);
            x = [A, -B'; B, S]\[zeros(fe.nrDoF,1); b];
        end

    end

end
classdef exBiharmonicHHJPatchTest < PatchTestCase
    %exBiharmonicHHJPatchTest - Patch test of the example exBiharmonicHHJ
    %   The Hellan-Herrmann-Johnson method HHJk x Pk(k+1) for the clamped
    %   plate reproduces the solution u = x^2(1-x)^2 y^2(1-y)^2 (degree 8)
    %   of problemBiharmonicSquareClampedExact and its Hessian for k = 7.
    %   Clamped polynomial solutions on the square have at least degree 8,
    %   so no lower degree gives a patch test for this example.
    %
    %   See also exBiharmonicHHJ, HHJk, problemBiharmonicSquareClampedExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        x
        fe
        pk
        p
    end

    methods (TestClassSetup)

        function solveProblem(testCase)
            testCase.p = problemBiharmonicSquareClampedExact();
            mesh = PatchTestCase.patchMesh(struct('name','Square'));
            [testCase.x,testCase.fe,testCase.pk] = exBiharmonicHHJPatchTest.solve(mesh,testCase.p,7);
        end

    end

    methods (Test)

        function testL2Error(testCase)
            p = testCase.p;
            u = testCase.x(testCase.fe.nrDoF+1:end);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(u,testCase.pk,16,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testH1Error(testCase)
            p = testCase.p;
            u = testCase.x(testCase.fe.nrDoF+1:end);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(u,testCase.pk,16,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testMomentError(testCase)
            p = testCase.p;
            sigma = testCase.x(1:testCase.fe.nrDoF);

            err = ErrorForm('( exact(x) - sigma )^2 dx','sigma').assemble(sigma,testCase.fe,16, ...
                'exact',@p.exactD2);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe,pk] = solve(mesh,p,degree)
            fe = HHJk(degree,mesh);
            pk = Pk(degree+1,mesh);
            A = BilinearForm('sigma:tau dx',{'sigma','tau'}).assemble(fe);
            B = BilinearForm('div(sigma)*grad(v) dx - [sigma*nt]*{grad(v)*t} ds',{'sigma','v'}).assemble(fe,pk);
            b = LinearForm('f(x)*v dx').assemble(pk,12,'f',@p.f);
            S = [A,B';B,sparse(pk.nrDoF,pk.nrDoF)];
            rhs = [zeros(fe.nrDoF,1);-b];
            freeNodes = [1:fe.nrDoF,fe.nrDoF+pk.freeNodes];
            x = zeros(size(S,1),1);
            x(freeNodes) = S(freeNodes,freeNodes)\rhs(freeNodes);
        end

    end

end
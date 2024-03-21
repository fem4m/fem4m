classdef exBiharmonicArgyrisPatchTest < PatchTestCase
    %exBiharmonicArgyrisPatchTest - Patch test of the example exBiharmonicArgyris
    %   The C1 conforming discretization of the clamped plate with the
    %   TUBA element TUBAk(8) (the alternative in the example) reproduces
    %   the solution u = x^2(1-x)^2 y^2(1-y)^2 (degree 8) of
    %   problemBiharmonicSquareClampedExact. The quintic Argyris element
    %   cannot represent any nonzero clamped polynomial on the square, so
    %   for Argyris only the C1 conformity of the discrete solution is
    %   checked (vanishing jumps of u and grad(u)).
    %
    %   See also exBiharmonicArgyris, TUBAk, Argyris,
    %   problemBiharmonicSquareClampedExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh
    end

    methods (TestClassSetup)

        function createMesh(testCase)
            testCase.mesh = PatchTestCase.patchMesh(struct('name','Square'));
        end

    end

    methods (Test)

        function testTUBAL2Error(testCase)
            p = problemBiharmonicSquareClampedExact();
            [x,fe] = exBiharmonicArgyrisPatchTest.solve(TUBAk(8,testCase.mesh),p);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,16,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testTUBAH1Error(testCase)
            p = problemBiharmonicSquareClampedExact();
            [x,fe] = exBiharmonicArgyrisPatchTest.solve(TUBAk(8,testCase.mesh),p);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,16,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testTUBAH2Error(testCase)
            p = problemBiharmonicSquareClampedExact();
            [x,fe] = exBiharmonicArgyrisPatchTest.solve(TUBAk(8,testCase.mesh),p);

            err = ErrorForm('( exact(x) - D2(u) )^2 dx').assemble(x,fe,16,'exact',@p.exactD2);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testArgyrisC1Conformity(testCase)
            p = problemBiharmonicSquareClampedExact();
            [x,fe] = exBiharmonicArgyrisPatchTest.solve(Argyris(5,testCase.mesh),p);

            eta = EstimatorForm('h^-1*([u])^2dsi + h^1*([grad(u)*n])^2dsi');
            eta4e = eta.assemble(x,fe);

            testCase.verifyLessThan(norm(eta4e),testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(fe,p)
            A = BilinearForm('D2(u):D2(v) dx').assemble(fe);
            b = LinearForm('f(x) v dx').assemble(fe,'f',@p.f);
            x = zeros(size(A,1),1);
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
        end

    end

end
classdef exPoissonP1PatchTest < PatchTestCase
    %exPoissonP1PatchTest - Patch test of the example exPoissonP1
    %   Conforming Pk finite elements for -Delta u = f with homogeneous
    %   Dirichlet data reproduce the biquadratic solution of
    %   problemScalarEllipticSquareExact for k >= 4.
    %
    %   See also exPoissonP1, problemScalarEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degree = struct('P4',4,'P5',5)
    end

    methods (Test)

        function testH1Error(testCase,degree)
            p = problemScalarEllipticSquareExact();
            [x,fe] = exPoissonP1PatchTest.solve(p,degree);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testL2Error(testCase,degree)
            p = problemScalarEllipticSquareExact();
            [x,fe] = exPoissonP1PatchTest.solve(p,degree);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(p,degree)
            mesh = PatchTestCase.patchMesh(struct('name','Square'));
            fe = Pk(degree,mesh);
            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);
            b = LinearForm('f(x) v dx').assemble(fe,'f',@p.f);
            x = zeros(size(A,1),1);
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
        end

    end

end
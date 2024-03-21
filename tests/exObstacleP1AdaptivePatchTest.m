classdef exObstacleP1AdaptivePatchTest < PatchTestCase
    %exObstacleP1AdaptivePatchTest - Patch test of the example exObstacleP1Adaptive
    %   The obstacle problem with Pk finite elements solved by the primal-
    %   dual active set method pdaSet reproduces the monomial
    %   u = x^i y^j, i+j <= k <= 3, for an inactive obstacle one unit
    %   below u with zero Lagrange multiplier. Higher degrees are not
    %   tested since pdaSet constrains only the nodal values.
    %
    %   See also exObstacleP1Adaptive, pdaSet,
    %   problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        monom = struct( ...
            'P1_x1',struct('degree',1,'i',1,'j',0), ...
            'P2_x1y1',struct('degree',2,'i',1,'j',1), ...
            'P3_x2y1',struct('degree',3,'i',2,'j',1))
    end

    methods (Test)

        function testInactiveObstacleError(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            mesh = PatchTestCase.patchMesh(struct('name','Square'));
            lower = @(x)(p.exactU(x)-1);
            [x,~,fe] = exObstacleP1AdaptivePatchTest.solve(mesh,p,monom.degree,lower,@(x)(inf(size(x,1),1)));

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testInactiveObstacleMultiplier(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            mesh = PatchTestCase.patchMesh(struct('name','Square'));
            lower = @(x)(p.exactU(x)-1);
            [~,lambda] = exObstacleP1AdaptivePatchTest.solve(mesh,p,monom.degree,lower,@(x)(inf(size(x,1),1)));

            testCase.verifyEqual(lambda,zeros(size(lambda)),'AbsTol',testCase.tol);
        end

    end

    methods (Static)

        function [x,lambda,fe] = solve(mesh,p,degree,lower,upper)
            fe = Pk(degree,mesh);
            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);
            b = LinearForm('f(x) v dx').assemble(fe,'f',@p.f);
            x = zeros(size(A,1),1);
            lambda = zeros(size(A,1),1);
            x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
            b = b - A*x;
            low = lower(fe.c4dof(fe.freeNodes,:));
            up = upper(fe.c4dof(fe.freeNodes,:));
            [x(fe.freeNodes),lambda(fe.freeNodes)] = pdaSet(A(fe.freeNodes,fe.freeNodes),b(fe.freeNodes),low,up);
        end

    end

end
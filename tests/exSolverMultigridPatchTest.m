classdef exSolverMultigridPatchTest < PatchTestCase
    %exSolverMultigridPatchTest - Patch test of the example exSolverMultigrid
    %   The multigrid iteration with Jacobi or Richardson smoothing for Pk
    %   finite elements on a hierarchy of uniformly refined meshes
    %   reproduces every monomial u = x^i y^j with i+j <= k (with
    %   inhomogeneous Dirichlet data) up to the solver tolerance.
    %
    %   See also exSolverMultigrid, Multigrid,
    %   problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (Constant)
        solverTol = 1e-10
        nrLevels = 4
    end

    properties (TestParameter)
        monom = struct( ...
            'P1_x1',struct('degree',1,'i',1,'j',0), ...
            'P1_y1',struct('degree',1,'i',0,'j',1), ...
            'P2_x1y1',struct('degree',2,'i',1,'j',1), ...
            'P3_x2y1',struct('degree',3,'i',2,'j',1))
        smoother = struct('Jacobi','Jacobi','Richardson','Richardson')
    end

    methods (Test)

        function testH1Error(testCase,monom,smoother)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,pk] = exSolverMultigridPatchTest.solve(p,monom.degree,smoother);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,pk,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,1e2*exSolverMultigridPatchTest.solverTol);
        end

    end

    methods (Static)

        function [x,pk] = solve(p,degree,smoother)
            mesh = Mesh(struct('name','Square'));
            mg = Multigrid(smoother);
            pk = [];
            for k = 1 : exSolverMultigridPatchTest.nrLevels
                mesh = mesh.refine('uniform');
                pkOld = pk;
                pk = Pk(degree,mesh);
                A = BilinearForm('grad(u) * grad(v) dx').assemble(pk);
                if k == 1
                    mg = mg.addLevel(A(pk.freeNodes,pk.freeNodes));
                else
                    P = pk.prolongation(pkOld);
                    mg = mg.addLevel(A(pk.freeNodes,pk.freeNodes),P(pk.freeNodes,pkOld.freeNodes));
                end
            end
            b = LinearForm('f(x) v dx').assemble(pk,'f',@p.f);
            x = zeros(size(A,1),1);
            x(pk.fixedNodes) = p.u4Db(pk.c4dof(pk.fixedNodes,:));
            b = b - A*x;
            x(pk.freeNodes) = mg.solve(b(pk.freeNodes),exSolverMultigridPatchTest.solverTol,100,[],0);
        end

    end

end
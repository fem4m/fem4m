classdef exSolverSchwarzPatchTest < PatchTestCase
    %exSolverSchwarzPatchTest - Patch test of the example exSolverSchwarz
    %   The preconditioned CG method with additive or multiplicative
    %   Schwarz preconditioners (HB, BPX, VO, BPS, NN, BDDC) for Pk finite
    %   elements on a hierarchy of uniformly refined meshes converges and
    %   reproduces the monomial u = x^i y^j with i+j <= k (with
    %   inhomogeneous Dirichlet data) up to the solver tolerance. The BPS
    %   preconditioner, which treats the edge, face and interior degrees
    %   of freedom separately, is also tested for P3 in 2d and in 3d on
    %   uniformly refined (bisection) tetrahedral meshes, where it
    %   reproduces a harmonic quadratic. The HB preconditioner marks as
    %   new dofs of a level exactly the free dofs at nodes that are no
    %   nodes of the previous level, also for P3 and P4, whose dof
    %   numbers change under refinement, and converges for them (P4 only
    %   additive: the default damping is too large for applyMult). On a
    %   surface mesh, whose sides are edges like in 2d, all
    %   preconditioners converge to the direct solution. For the
    %   symmetric interior penalty method with the
    %   discontinuous elements DGDk and DGPk in 2d and 3d, VO, whose
    %   patches then hold all dofs of the patch elements, converges to the
    %   direct solution.
    %
    %   See also exSolverSchwarz, Schwarz,
    %   problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (Constant)
        solverTol = 1e-10
        nrLevels = 4
    end

    properties (TestParameter)
        preconditioner = struct('HB','SchwarzHB','BPX','SchwarzBPX','VO','SchwarzVO', ...
            'BPS','SchwarzBPS','NN','SchwarzNN','BDDC','SchwarzBDDC')
        type = struct('add','add','mult','mult')
        monom = struct( ...
            'P1_x1',struct('degree',1,'i',1,'j',0), ...
            'P2_x1y1',struct('degree',2,'i',1,'j',1))
        hbDegree = struct('P1',1,'P2',2,'P3',3,'P4',4)
        % applyMult needs omega*lambda_max(A) < 2, which the default
        % damping 1/8 of SchwarzHB violates from degree 4 on
        hbMonom = struct( ...
            'P3_x2y1_add',struct('degree',3,'i',2,'j',1,'type','add'), ...
            'P3_x2y1_mult',struct('degree',3,'i',2,'j',1,'type','mult'), ...
            'P4_x3y1_add',struct('degree',4,'i',3,'j',1,'type','add'))
        discontinuous = struct( ...
            'DGDk',struct('name','DGDk','degree',1,'dim',2), ...
            'DGPk',struct('name','DGPk','degree',2,'dim',2), ...
            'DGDk3d',struct('name','DGDk3d','degree',1,'dim',3), ...
            'DGPk3d',struct('name','DGPk3d','degree',1,'dim',3))
    end

    methods (Test)

        function testConvergence(testCase,preconditioner,type,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [~,~,flag] = exSolverSchwarzPatchTest.solve(p,monom.degree,preconditioner,type);

            testCase.verifyEqual(flag,0);
        end

        function testH1Error(testCase,preconditioner,type,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,pk] = exSolverSchwarzPatchTest.solve(p,monom.degree,preconditioner,type);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,pk,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,1e2*exSolverSchwarzPatchTest.solverTol);
        end

        function testBPSHigherOrder2d(testCase,type)
            p = problemScalarEllipticSquareMonomExact(2,1);
            [x,pk,flag] = exSolverSchwarzPatchTest.solve(p,3,'SchwarzBPS',type);

            testCase.verifyEqual(flag,0);
            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,pk,'exact',@p.exactGradU);
            testCase.verifyLessThan(err,1e2*exSolverSchwarzPatchTest.solverTol);
        end

        function testHBNewDofs(testCase,hbDegree)
            % the new dofs are those at the nodes of the fine mesh that
            % are no nodes of the coarse mesh, whatever their numbers
            mesh = Mesh(struct('name','Square'));
            mesh = mesh.refine('uniform');
            pkOld = Pk(hbDegree,mesh);
            a = BilinearForm('grad(u) * grad(v) dx');
            schwarz = SchwarzHB().addLevel(pkOld,[],[],a.assemble(pkOld));
            mesh = mesh.refine('uniform');
            pk = Pk(hbDegree,mesh);
            schwarz = schwarz.addLevel(pk,pkOld,pk.prolongation(pkOld),a.assemble(pk));

            atOldNode = ismembertol(pk.c4dof,pkOld.c4dof,1e-12,'ByRows',true);
            expected = schwarz.omega*double(~atOldNode(pk.freeNodes));
            testCase.verifyEqual(full(diag(schwarz.B{end})),expected);
        end

        function testHBHigherOrder(testCase,hbMonom)
            p = problemScalarEllipticSquareMonomExact(hbMonom.i,hbMonom.j);
            [x,pk,flag] = exSolverSchwarzPatchTest.solve(p,hbMonom.degree,'SchwarzHB',hbMonom.type);

            testCase.verifyEqual(flag,0);
            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,pk,'exact',@p.exactGradU);
            testCase.verifyLessThan(err,1e2*exSolverSchwarzPatchTest.solverTol);
        end

        function testBPS3d(testCase,type)
            % the harmonic quadratic u = x*y + z with f = 0
            u = @(z)z(:,1).*z(:,2) + z(:,3);
            mesh = Mesh3d(Geometry(struct('name','Cube')));
            schwarz = SchwarzBPS(0.4);
            pk = [];
            for k = 1 : 3
                if k > 1
                    mesh = mesh.refine('uniform');
                end
                pkOld = pk;
                pk = Pk3d(3,mesh);
                [A,A4e] = BilinearForm('grad(u) * grad(v) dx').assemble(pk);
                if k == 1
                    P = [];
                else
                    P = pk.prolongation(pkOld);
                end
                schwarz = schwarz.addLevel(pk,pkOld,P,A,A4e);
            end
            x = pk.project(u);
            expected = x;
            x(pk.freeNodes) = 0;
            b = -A*x;
            if strcmp(type,'add')
                B = @(b)(schwarz.applyAdd(b));
            else
                B = @(b)(schwarz.applyMult(b));
            end

            [x(pk.freeNodes),flag] = pcg(A(pk.freeNodes,pk.freeNodes),b(pk.freeNodes), ...
                exSolverSchwarzPatchTest.solverTol,1000,B);

            testCase.verifyEqual(flag,0);
            testCase.verifyEqual(x,expected,'AbsTol',1e3*exSolverSchwarzPatchTest.solverTol);
        end

        function testSurface(testCase,preconditioner)
            p = problemSurfaceSquare();
            mesh = MeshSurface(p.geom,@p.surfaceFunction,@p.gradSurfaceFunction,[]);
            schwarz = feval(preconditioner);
            pk = [];
            for k = 1 : 3
                mesh = mesh.refine('uniform');
                pkOld = pk;
                pk = Pk(2,mesh);
                [A,A4e] = BilinearForm('grad(u) * grad(v) dx').assemble(pk);
                if k == 1
                    P = [];
                else
                    P = pk.prolongation(pkOld);
                end
                schwarz = schwarz.addLevel(pk,pkOld,P,A,A4e);
            end
            x = pk.project(@(z)z(:,1).*z(:,2) + z(:,3));
            x(pk.freeNodes) = 0;
            b = -A*x;
            expected = x;
            expected(pk.freeNodes) = A(pk.freeNodes,pk.freeNodes)\b(pk.freeNodes);

            [x(pk.freeNodes),flag] = pcg(A(pk.freeNodes,pk.freeNodes),b(pk.freeNodes), ...
                exSolverSchwarzPatchTest.solverTol,1000,@(r)schwarz.applyAdd(r));

            testCase.verifyEqual(flag,0);
            testCase.verifyLessThan(norm(x-expected)/norm(expected), ...
                1e3*exSolverSchwarzPatchTest.solverTol);
        end

        function testDiscontinuousVO(testCase,discontinuous)
            % the patches of the discontinuous elements hold all dofs of
            % the patch elements; symmetric interior penalty system
            if discontinuous.dim == 2
                mesh = Mesh(struct('name','Square'));
            else
                mesh = Mesh3d(Geometry(struct('name','Cube')));
            end
            degree = discontinuous.degree;
            schwarz = SchwarzVO();
            fe = [];
            for k = 1 : 2
                mesh = mesh.refine('uniform');
                feOld = fe;
                fe = feval(discontinuous.name,degree,mesh);
                A = BilinearForm('grad(u)*grad(v)dx').assemble(fe);
                C = BilinearForm('[u]{grad(v).n}ds').assemble(fe);
                D = BilinearForm('h^-1*penalty(x)*[u][v] ds').assemble(fe,[],'penalty',10*degree^2);
                A = A - C - C' + D;
                if k == 1
                    P = [];
                else
                    P = fe.prolongation(feOld);
                end
                schwarz = schwarz.addLevel(fe,feOld,P,A);
            end
            b = LinearForm('f(x)*v dx').assemble(fe,'f',@(x)1+x(:,1));
            expected = A\b;

            [x,flag] = pcg(A,b,exSolverSchwarzPatchTest.solverTol,1000,@(r)schwarz.applyAdd(r));

            testCase.verifyEqual(flag,0);
            testCase.verifyLessThan(norm(x-expected)/norm(expected), ...
                1e3*exSolverSchwarzPatchTest.solverTol);
        end


    end

    methods (Static)

        function [x,pk,flag] = solve(p,degree,preconditioner,type)
            mesh = Mesh(struct('name','Square'));
            schwarz = feval(preconditioner);
            pk = [];
            for k = 1 : exSolverSchwarzPatchTest.nrLevels
                mesh = mesh.refine('uniform');
                pkOld = pk;
                pk = Pk(degree,mesh);
                [A,A4e] = BilinearForm('grad(u) * grad(v) dx').assemble(pk);
                if k == 1
                    P = [];
                else
                    P = pk.prolongation(pkOld);
                end
                schwarz = schwarz.addLevel(pk,pkOld,P,A,A4e);
            end
            if strcmpi(type,'add')
                B = @(b)(schwarz.applyAdd(b));
            else
                B = @(b)(schwarz.applyMult(b));
            end
            b = LinearForm('f(x) v dx').assemble(pk,'f',@p.f);
            x = zeros(size(A,1),1);
            x(pk.fixedNodes) = p.u4Db(pk.c4dof(pk.fixedNodes,:));
            b = b - A*x;
            [x(pk.freeNodes),flag] = pcg(A(pk.freeNodes,pk.freeNodes),b(pk.freeNodes), ...
                exSolverSchwarzPatchTest.solverTol,1000,B);
        end

    end

end
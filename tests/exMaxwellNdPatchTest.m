classdef exMaxwellNdPatchTest < PatchTestCase
    %exMaxwellNdPatchTest - Patch test of the example exMaxwellNd
    %   Nedelec elements for curl(rot(u)) + alpha*u = f with inhomogeneous
    %   tangential Dirichlet data reproduce the field
    %   u = (x^i y^j, x^j y^i) and its rot for i+j <= k-1 with Nd1k(k)
    %   (first kind) and for i+j <= k with Nd2k(k) (second kind).
    %
    %   See also exMaxwellNd, Nd1k, Nd2k, problemMaxwellSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (ClassSetupParameter)
        domain = struct('Square',struct('name','Square'), ...
            'Lshape',struct('name','Lshape'))
    end

    properties (TestParameter)
        monom = struct( ...
            'Nd1k1_const',struct('space','Nd1k','degree',1,'i',0,'j',0), ...
            'Nd1k2_x1',struct('space','Nd1k','degree',2,'i',1,'j',0), ...
            'Nd1k3_x1y1',struct('space','Nd1k','degree',3,'i',1,'j',1), ...
            'Nd1k3_x2',struct('space','Nd1k','degree',3,'i',2,'j',0), ...
            'Nd2k1_x1',struct('space','Nd2k','degree',1,'i',1,'j',0), ...
            'Nd2k2_x1y1',struct('space','Nd2k','degree',2,'i',1,'j',1), ...
            'Nd2k2_y2',struct('space','Nd2k','degree',2,'i',0,'j',2), ...
            'Nd2k3_x2y1',struct('space','Nd2k','degree',3,'i',2,'j',1))
    end

    properties
        mesh
    end

    methods (TestClassSetup)

        function createMesh(testCase,domain)
            testCase.mesh = PatchTestCase.patchMesh(domain);
        end

    end

    methods (Test)

        function testL2Error(testCase,monom)
            p = problemMaxwellSquareMonomExact(1,monom.i,monom.j);
            [x,nd] = exMaxwellNdPatchTest.solve(testCase.mesh,p,monom);

            err = ErrorForm('(exact(x)-u)^2 dx').assemble(x,nd,12,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testRotError(testCase,monom)
            p = problemMaxwellSquareMonomExact(1,monom.i,monom.j);
            [x,nd] = exMaxwellNdPatchTest.solve(testCase.mesh,p,monom);

            err = ErrorForm('(exact(x)-rot(u))^2 dx').assemble(x,nd,12,'exact',@p.exactRotU);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,nd] = solve(mesh,p,monom)
            nd = feval(monom.space,monom.degree,mesh);
            A = BilinearForm('rot(u)*rot(v) dx + alpha(x)*u*v dx').assemble(nd,'alpha',p.alpha);
            b = LinearForm('f(x)*v dx').assemble(nd,'f',@p.f);
            x = nd.projectBoundary(mesh.ed4Db,@(pts,tangents)sum(p.u4Db(pts).*tangents,2));
            b = b - A*x;
            free = nd.freeNodes;
            x(free) = A(free,free)\b(free);
        end

    end

end
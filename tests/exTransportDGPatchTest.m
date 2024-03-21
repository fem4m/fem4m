classdef exTransportDGPatchTest < PatchTestCase
    %exTransportDGPatchTest - Patch test of the example exTransportDG
    %   The dG method DGPk for the transport equation beta*grad(u) = f with
    %   centered (upwind = 0) or upwind (upwind = 1) fluxes and inflow data
    %   imposed weakly reproduces every monomial u = x^i y^j with
    %   i+j <= k, including nonzero inflow data. With k = 4 it also
    %   reproduces the solution of problemTransportSquareExact.
    %
    %   See also exTransportDG, problemTransportSquareMonomExact,
    %   problemTransportSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        monom = struct( ...
            'DGP1_x1',struct('degree',1,'i',1,'j',0), ...
            'DGP1_y1',struct('degree',1,'i',0,'j',1), ...
            'DGP2_x1y1',struct('degree',2,'i',1,'j',1), ...
            'DGP2_x2',struct('degree',2,'i',2,'j',0), ...
            'DGP3_x2y1',struct('degree',3,'i',2,'j',1), ...
            'DGP3_y3',struct('degree',3,'i',0,'j',3))
        upwind = struct('centered',0,'upwind',1)
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

        function testL2Error(testCase,monom,upwind)
            p = problemTransportSquareMonomExact(monom.i,monom.j);
            [x,fe] = exTransportDGPatchTest.solve(testCase.mesh,p,monom.degree,upwind);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,12,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testSquareExactP4(testCase,upwind)
            p = problemTransportSquareExact();
            [x,fe] = exTransportDGPatchTest.solve(testCase.mesh,p,4,upwind);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,12,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(mesh,p,degree,upwind)
            normal4ed = mesh.normal4ed;
            betan = @(pts,parts)(sum(p.beta(pts).*normal4ed(parts,:),2));
            fe = DGPk(degree,mesh);
            a = BilinearForm(['beta(x)*grad(u)vdx - betan(x)*[u]{v} dsi', ...
                ' + betanm(x)*[u]*{v} dsd + delta(x)*[u][v]dsi']);
            A = a.assemble(fe,'beta',@p.beta, ...
                'betan',@(pts,pts_ref,parts)(betan(pts,parts)), ...
                'betanm',@(pts,pts_ref,parts)(0.5*(abs(betan(pts,parts))-betan(pts,parts))), ...
                'delta',@(pts,pts_ref,parts)(upwind/2*abs(betan(pts,parts))));
            b = LinearForm('f(x)*v dx + betanm(x)*v dsd');
            b = b.assemble(fe,'f',@p.f,'betanm',@(pts,pts_ref,parts,normals) ...
                (0.5*p.u4Db(pts).*(abs(sum(p.beta(pts).*normals,2))-sum(p.beta(pts).*normals,2))));
            x = A\b;
        end

    end

end
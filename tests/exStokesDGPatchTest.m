classdef exStokesDGPatchTest < PatchTestCase
    %exStokesDGPatchTest - Patch test of the example exStokesDG
    %   The equal-order symmetric interior penalty dG method
    %   DGPkVector x DGPk with pressure jump stabilization and
    %   inhomogeneous Dirichlet data imposed by Nitsche's method reproduces
    %   the divergence-free velocity u = curl(x^i y^j) in P_k and the
    %   continuous mean-zero pressure p = x^m - 1/(m+1) with m <= k.
    %
    %   See also exStokesDG, problemStokesSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        flow = struct( ...
            'DGP1_psix1y1_p0',struct('degree',1,'i',1,'j',1,'m',0), ...
            'DGP1_psix2_px1',struct('degree',1,'i',2,'j',0,'m',1), ...
            'DGP2_psix2y1_px2',struct('degree',2,'i',2,'j',1,'m',2), ...
            'DGP2_psix1y2_px1',struct('degree',2,'i',1,'j',2,'m',1), ...
            'DGP2_psix3_px2',struct('degree',2,'i',3,'j',0,'m',2), ...
            'DGP3_psix2y2_px3',struct('degree',3,'i',2,'j',2,'m',3))
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

        function testVelocityError(testCase,flow)
            p = problemStokesSquareMonomExact(1,flow.i,flow.j,flow.m);
            [x,feV] = exStokesDGPatchTest.solve(testCase.mesh,p,flow.degree);

            err = ErrorForm('(exact(x)-u)^2 dx').assemble(x(1:feV.nrDoF),feV,6,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testPressureError(testCase,flow)
            p = problemStokesSquareMonomExact(1,flow.i,flow.j,flow.m);
            [x,feV,feP] = exStokesDGPatchTest.solve(testCase.mesh,p,flow.degree);

            err = ErrorForm('(exact(x)-u)^2 dx').assemble(x(feV.nrDoF+1:end),feP,6,'exact',@p.exactP);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,feV,feP] = solve(mesh,p,degree)
            penalty = 2*(degree+1)*degree;
            feV = DGPkVector(degree,mesh);
            feP = DGPk(degree,mesh);
            a = BilinearForm(['grad(u):grad(v)dx', ...
                '- [uxn]:{grad(v)}ds - {grad(u)}:[vxn]ds', ...
                '+ h^-1*alpha(x)*[uxn]:[vxn] ds']);
            A = a.assemble(feV,'alpha',penalty);
            B = BilinearForm(['-div(v)*q dx','+ [v*n]{q} ds'],{'v','q'}).assemble(feV,feP);
            C = LinearForm('v dx').assemble(feP);
            S = BilinearForm('h^1*[p][q] dsi',{'p','q'}).assemble(feP);
            O = sparse(1,size(A,2));
            b1 = LinearForm('f(x)*v dx').assemble(feV,'f',@p.f);
            rhsNitsche = LinearForm('h^-1*p(x)*v dsd - d(x)*grad(v)n dsd');
            b1 = b1 + rhsNitsche.assemble(feV,2*degree+2,'p',@(x)(penalty*p.u4Db(x)),'d',@(x)(p.u4Db(x)));
            b2 = LinearForm('d(x)*vn dsd').assemble(feP,2*degree+2,'d',@(x)(p.u4Db(x)));
            x = [A,B',O'; -B,S,C; O,C',0]\[b1;-b2;0];
            x = x(1:end-1);
        end

    end

end
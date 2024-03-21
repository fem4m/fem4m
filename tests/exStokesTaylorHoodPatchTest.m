classdef exStokesTaylorHoodPatchTest < PatchTestCase
    %exStokesTaylorHoodPatchTest - Patch test of the example exStokesTaylorHood
    %   Conforming Stokes discretizations with inhomogeneous Dirichlet data
    %   and a mean-zero pressure constraint reproduce the divergence-free
    %   velocity u = curl(x^i y^j) and the pressure p = x^m - 1/(m+1)
    %   whenever both lie in the discrete spaces:
    %       Taylor-Hood P2-P1   u in P2, p in P1
    %       Taylor-Hood P3-P2   u in P3, p in P2
    %       MINI P1+bubble-P1   u in P1, p in P1
    %       Crouzeix-Raviart-P0 u in P1, p = 0
    %
    %   See also exStokesTaylorHood, problemStokesSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (ClassSetupParameter)
        reftype = struct('criss','criss','crisscross','crisscross')
    end

    properties (TestParameter)
        flow = struct( ...
            'TH2_psix1y1_p0',struct('element','TH2','i',1,'j',1,'m',0), ...
            'TH2_psix2y1_px1',struct('element','TH2','i',2,'j',1,'m',1), ...
            'TH2_psix3_px1',struct('element','TH2','i',3,'j',0,'m',1), ...
            'TH3_psix3y1_px2',struct('element','TH3','i',3,'j',1,'m',2), ...
            'TH3_psix2y2_px1',struct('element','TH3','i',2,'j',2,'m',1), ...
            'TH3_psix4_px2',struct('element','TH3','i',4,'j',0,'m',2), ...
            'MINI_psix1y1_px1',struct('element','MINI','i',1,'j',1,'m',1), ...
            'MINI_psix2_px1',struct('element','MINI','i',2,'j',0,'m',1), ...
            'CR_psix1y1_p0',struct('element','CR','i',1,'j',1,'m',0), ...
            'CR_psiy2_p0',struct('element','CR','i',0,'j',2,'m',0))
    end

    properties
        mesh
    end

    methods (TestClassSetup)

        function createMesh(testCase,reftype)
            testCase.mesh = PatchTestCase.patchMesh(struct('name','Square','reftype',reftype));
        end

    end

    methods (Test)

        function testVelocityError(testCase,flow)
            p = problemStokesSquareMonomExact(1,flow.i,flow.j,flow.m);
            [x,feV] = exStokesTaylorHoodPatchTest.solve(testCase.mesh,p,flow.element);

            err = ErrorForm('(exact(x)-u)^2 dx').assemble(x(1:feV.nrDoF),feV,6,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testPressureError(testCase,flow)
            p = problemStokesSquareMonomExact(1,flow.i,flow.j,flow.m);
            [x,feV,feP] = exStokesTaylorHoodPatchTest.solve(testCase.mesh,p,flow.element);

            err = ErrorForm('(exact(x)-u)^2 dx').assemble(x(feV.nrDoF+1:end),feP,6,'exact',@p.exactP);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [feV,feP] = spaces(mesh,element)
            switch element
                case 'TH2'
                    feV = PkVector(2,mesh);
                    feP = Pk(1,mesh);
                case 'TH3'
                    feV = PkVector(3,mesh);
                    feP = Pk(2,mesh);
                case 'MINI'
                    feV = PkVector(1,mesh) + FEVector(Bk(3,mesh),Bk(3,mesh));
                    feP = Pk(1,mesh);
                case 'CR'
                    feV = FEVector(CRk(1,mesh),CRk(1,mesh));
                    feP = DGDk(0,mesh);
            end
        end

        function [x,feV,feP] = solve(mesh,p,element)
            [feV,feP] = exStokesTaylorHoodPatchTest.spaces(mesh,element);
            A = BilinearForm('grad(u)*grad(v) dx').assemble(feV);
            B = BilinearForm('-div(u)*v dx').assemble(feV,feP);
            c = LinearForm('v dx').assemble(feP);
            b = LinearForm('f(x)*v dx').assemble(feV,'f',@p.f);
            x = zeros(feV.nrDoF+feP.nrDoF+1,1);
            if isprop(feV.fe{1},'c4dof')
                x(feV.fixedNodes) = reshape(p.u4Db(feV.fe{1}.c4dof(feV.fe{1}.fixedNodes,:)),[],1);
            else
                mid4ed = mesh.point4ed(1/2);
                x(feV.fixedNodes) = reshape(p.u4Db(mid4ed(feV.fe{1}.fixedNodes,:)),[],1);
            end
            S = [A,B',sparse(feV.nrDoF,1);B,sparse(feP.nrDoF,feP.nrDoF),c;sparse(1,feV.nrDoF),c',0];
            rhs = [b;zeros(feP.nrDoF+1,1)] - S*x;
            free = [feV.freeNodes(:);feV.nrDoF+(1:feP.nrDoF+1)'];
            x(free) = S(free,free)\rhs(free);
            x = x(1:end-1);
        end

    end

end
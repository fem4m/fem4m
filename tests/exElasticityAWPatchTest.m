classdef exElasticityAWPatchTest < PatchTestCase
    %exElasticityAWPatchTest - Patch test of the example exElasticityAW
    %   The Arnold-Winther mixed method AWk(k+2) x DGDkVector(k) with natural
    %   displacement and essential traction data reproduces the
    %   displacement u = (x^i y^j, x^j y^i) and its stress for i+j <= k.
    %   The stress is also reproduced for i+j = k+1 and i+j = k+2 (sigma
    %   in P_{k+1}, div(sigma) in P_k), while u_h is then only an
    %   approximation, with pure displacement data and with traction
    %   data, whose vertex constraints (Constraints and the second output
    %   bN of interpolateBoundary) the example imposes with a multiplier.
    %
    %   See also exElasticityAW, AWk, problemElasticitySquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (ClassSetupParameter)
        domain = struct('Square',struct('name','Square','bdry','DNDN'), ...
            'Lshape',struct('name','Lshape','bdry','DN'))
    end

    properties (TestParameter)
        monom = struct( ...
            'AW1_x1',struct('degree',1,'i',1,'j',0), ...
            'AW1_y1',struct('degree',1,'i',0,'j',1), ...
            'AW2_x2',struct('degree',2,'i',2,'j',0), ...
            'AW2_x1y1',struct('degree',2,'i',1,'j',1))
        stressMonom = struct( ...
            'AW1_x1y1',struct('degree',1,'i',1,'j',1), ...
            'AW1_y2',struct('degree',1,'i',0,'j',2), ...
            'AW1_x2y1',struct('degree',1,'i',2,'j',1), ...
            'AW1_x3',struct('degree',1,'i',3,'j',0), ...
            'AW2_x2y1',struct('degree',2,'i',2,'j',1), ...
            'AW2_y3',struct('degree',2,'i',0,'j',3), ...
            'AW2_x3y1',struct('degree',2,'i',3,'j',1))
    end

    properties
        mesh
        dirichletMesh
    end

    methods (TestClassSetup)

        function createMesh(testCase,domain)
            testCase.mesh = PatchTestCase.patchMesh(domain);
            domain.bdry = 'D';
            testCase.dirichletMesh = PatchTestCase.patchMesh(domain);
        end

    end

    methods (Test)

        function testDisplacementError(testCase,monom)
            p = problemElasticitySquareMonomExact('E',2,'nu',0.3,monom.i,monom.j);
            [x,fe,dg] = exElasticityAWPatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x(fe.nrDoF+1:end),dg,12,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testStressError(testCase,monom)
            p = problemElasticitySquareMonomExact('E',2,'nu',0.3,monom.i,monom.j);
            [x,fe] = exElasticityAWPatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - sigma )^2 Cin(x) dx','sigma').assemble(x(1:fe.nrDoF),fe,12, ...
                'exact',@p.exactSigma,'Cin',p.CInv);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testStressErrorDirichlet(testCase,stressMonom)
            p = problemElasticitySquareMonomExact('E',2,'nu',0.3,stressMonom.i,stressMonom.j);
            [x,fe] = exElasticityAWPatchTest.solve(testCase.dirichletMesh,p,stressMonom.degree);

            err = ErrorForm('( exact(x) - sigma )^2 Cin(x) dx','sigma').assemble(x(1:fe.nrDoF),fe,12, ...
                'exact',@p.exactSigma,'Cin',p.CInv);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testStressErrorTraction(testCase,stressMonom)
            p = problemElasticitySquareMonomExact('E',2,'nu',0.3,stressMonom.i,stressMonom.j);
            [x,fe] = exElasticityAWPatchTest.solve(testCase.mesh,p,stressMonom.degree);

            err = ErrorForm('( exact(x) - sigma )^2 Cin(x) dx','sigma').assemble(x(1:fe.nrDoF),fe,12, ...
                'exact',@p.exactSigma,'Cin',p.CInv);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe,dg] = solve(mesh,p,degree)
            fe = AWk(degree+2,mesh);
            dg = DGDkVector(degree,mesh);
            A = BilinearForm('Cinv(x)*sigma*tau dx',{'sigma','tau'}).assemble(fe,'Cinv',p.CInv);
            B = BilinearForm('div(sigma)*v dx',{'sigma','v'}).assemble(fe,dg);
            bD = LinearForm('uD(x)*(tau*n) dsd',{'sigma','tau'}).assemble(fe,'uD',@p.u4Db);
            b = LinearForm('f(x)*v dx').assemble(dg,'f',@p.f);
            [x,bN] = fe.interpolateBoundary(@p.g);
            C = fe.Constraints;
            nrC = size(C,1);
            S = [A,B',C';B,sparse(dg.nrDoF,dg.nrDoF+nrC);C,sparse(nrC,dg.nrDoF+nrC)];
            x = [x;zeros(dg.nrDoF+nrC,1)];
            rhs = [bD;-b;bN] - S*x;
            free = [fe.freeNodes,fe.nrDoF+1:size(S,1)];
            x(free) = S(free,free)\rhs(free);
            x = x(1:fe.nrDoF+dg.nrDoF);
        end

    end

end
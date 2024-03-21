classdef exPoissonRTAdaptivePatchTest < PatchTestCase
    %exPoissonRTAdaptivePatchTest - Patch test of the example exPoissonRTAdaptive
    %   The Raviart-Thomas mixed method RTk x DGDk(k-1) with natural
    %   Dirichlet and essential Neumann data reproduces the flux
    %   sigma = grad(u) of every monomial u = x^i y^j with i+j <= k, and
    %   the discrete solution u_h is the L2 projection of u onto
    %   DGDk(k-1). In particular, u_h = u for i+j <= k-1. The error
    %   estimator vanishes.
    %
    %   See also exPoissonRTAdaptive, RTk,
    %   problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (ClassSetupParameter)
        domain = struct('Square',struct('name','Square','bdry','DNDN'), ...
            'Lshape',struct('name','Lshape','bdry','DN'))
    end

    properties (TestParameter)
        monom = struct( ...
            'RT1_x1',struct('degree',1,'i',1,'j',0), ...
            'RT2_x1y1',struct('degree',2,'i',1,'j',1), ...
            'RT2_y2',struct('degree',2,'i',0,'j',2), ...
            'RT3_x2y1',struct('degree',3,'i',2,'j',1), ...
            'RT3_y3',struct('degree',3,'i',0,'j',3))
        lowMonom = struct( ...
            'RT1_const',struct('degree',1,'i',0,'j',0), ...
            'RT2_y1',struct('degree',2,'i',0,'j',1), ...
            'RT3_x1y1',struct('degree',3,'i',1,'j',1), ...
            'RT3_x2',struct('degree',3,'i',2,'j',0))
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

        function testSigmaError(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,rt] = exPoissonRTAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('(exact(x)-sigma)^2 dx','sigma').assemble(x,rt, ...
                max(12,2*(monom.degree+1)),'exact',@p.exactSigma);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testUIsL2Projection(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,rt,dg] = exPoissonRTAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            uProjected = dg.project(@p.exactU);

            testCase.verifyEqual(x(rt.nrDoF+1:end),uProjected,'AbsTol',testCase.tol);
        end

        function testUError(testCase,lowMonom)
            p = problemScalarEllipticSquareMonomExact(lowMonom.i,lowMonom.j);
            [x,rt,dg] = exPoissonRTAdaptivePatchTest.solve(testCase.mesh,p,lowMonom.degree);

            err = ErrorForm('(exact(x)-u)^2 dx').assemble(x(rt.nrDoF+1:end),dg, ...
                max(12,2*(lowMonom.degree+1)),'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testEstimator(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,rt] = exPoissonRTAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            pkV = DGPkVector(monom.degree,testCase.mesh);
            xV = pkV.project(@(pts,ptsref,ind)(rt.getValue(pts,ptsref,ind,x(1:rt.nrDoF))));
            eta = EstimatorForm(['h^2*(f(x)+div(sigma))^2 dx + h^2*(rot(sigma))^2 dx', ...
                ' + h^1*([sigma*t])^2 dsi'],'sigma');
            eta4e = eta.assemble(xV,pkV,'f',@p.f);

            testCase.verifyLessThan(norm(eta4e),testCase.tol);
        end

    end

    methods (Static)

        function [x,rt,dg] = solve(mesh,p,degree)
            rt = RTk(degree,mesh);
            dg = DGDk(degree-1,mesh);
            A = BilinearForm('sigma*tau dx',{'sigma','tau'}).assemble(rt);
            B = BilinearForm('div(sigma)*v dx',{'sigma','v'}).assemble(rt,dg);
            b = LinearForm('f(x)*v dx').assemble(dg,'f',@(z)(-p.f(z)));
            uD = LinearForm('bd(x)*tau.n dsd',{'sigma','tau'}).assemble(rt,'bd',@p.u4Db);
            x = zeros(rt.nrDoF+dg.nrDoF,1);
            x(1:rt.nrDoF) = rt.projectNbBoundary(@p.g);
            b = [uD;b] - [A;B]*x(1:rt.nrDoF);
            freeNodes = [rt.freeNodes,rt.nrDoF+dg.freeNodes];
            S = [A(rt.freeNodes,rt.freeNodes),B(:,rt.freeNodes)';
                B(:,rt.freeNodes),sparse(dg.nrDoF,dg.nrDoF)];
            x(freeNodes) = S\b(freeNodes);
        end

    end

end
classdef exElasticityPkAdaptivePatchTest < PatchTestCase
    %exElasticityPkAdaptivePatchTest - Patch test of the example exElasticityPkAdaptive
    %   Conforming PkVector finite elements for linear elasticity with
    %   inhomogeneous displacement and traction data reproduce every
    %   displacement u = (x^i y^j, x^j y^i) with i+j <= k, and the
    %   residual-based error estimator vanishes. With k = 4 they also
    %   reproduce the solution of problemElasticitySquareExact.
    %
    %   See also exElasticityPkAdaptive, problemElasticitySquareMonomExact,
    %   problemElasticitySquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties (ClassSetupParameter)
        domain = struct('Square',struct('name','Square','bdry','DNDN'), ...
            'Lshape',struct('name','Lshape','bdry','DN'))
    end

    properties (TestParameter)
        monom = struct( ...
            'P1_x1',struct('degree',1,'i',1,'j',0), ...
            'P2_x1y1',struct('degree',2,'i',1,'j',1), ...
            'P2_y2',struct('degree',2,'i',0,'j',2), ...
            'P3_x2y1',struct('degree',3,'i',2,'j',1), ...
            'P3_y3',struct('degree',3,'i',0,'j',3))
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

        function testDisplacementError(testCase,monom)
            p = problemElasticitySquareMonomExact('E',2,'nu',0.3,monom.i,monom.j);
            [x,fe] = exElasticityPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('(exact(x)-u)^2 dx').assemble(x,fe,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testStressError(testCase,monom)
            p = problemElasticitySquareMonomExact('E',2,'nu',0.3,monom.i,monom.j);
            [x,fe] = exElasticityPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('(exact(x)-C(x)*eps(u))^2 Cin(x) dx').assemble(x,fe, ...
                'C',p.C,'exact',@p.exactSigma,'Cin',p.CInv);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testEstimator(testCase,monom)
            p = problemElasticitySquareMonomExact('E',2,'nu',0.3,monom.i,monom.j);
            [x,fe] = exElasticityPkAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            eta = EstimatorForm(['h^2*(f(x)+div(C(x)*eps(u)))^2 dx', ...
                '+ h^1*([C(x)*eps(u)*n])^2 dsi', ...
                '+ h^1*(g(x)-C(x)*eps(u)*n)^2 dsn']);
            eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n),'C',p.C);

            testCase.verifyLessThan(norm(eta4e),testCase.tol);
        end

        function testSquareExactP4(testCase)
            p = problemElasticitySquareExact('E',1e5,'nu',0.3);
            [x,fe] = exElasticityPkAdaptivePatchTest.solve(testCase.mesh,p,4);

            err = ErrorForm('(exact(x)-u)^2 dx').assemble(x,fe,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(mesh,p,degree)
            fe = PkVector(degree,mesh);
            A = BilinearForm('C(x)*eps(u):eps(v) dx').assemble(fe,'C',p.C);
            b = LinearForm('f(x)*v dx + g(x)*v dsn').assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
            x = zeros(fe.nrDoF,1);
            u4Db = fe.interpolate(@p.u4Db);
            x(fe.fixedNodes) = u4Db(fe.fixedNodes);
            b = b - A*x;
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
        end

    end

end
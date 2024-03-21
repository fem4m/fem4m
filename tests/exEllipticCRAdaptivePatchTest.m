classdef exEllipticCRAdaptivePatchTest < PatchTestCase
    %exEllipticCRAdaptivePatchTest - Patch test of the example exEllipticCRAdaptive
    %   Crouzeix-Raviart elements CRk (k = 1,3) with inhomogeneous Dirichlet
    %   and Neumann data reproduce every monomial u = x^i y^j with
    %   i+j <= k.
    %
    %   See also exEllipticCRAdaptive, CRk,
    %   problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties (ClassSetupParameter)
        domain = struct('Square',struct('name','Square','bdry','DNDN'), ...
            'Lshape',struct('name','Lshape','bdry','DN'))
    end

    properties (TestParameter)
        monom = struct( ...
            'CR1_x1',struct('degree',1,'i',1,'j',0), ...
            'CR1_y1',struct('degree',1,'i',0,'j',1), ...
            'CR3_x1y1',struct('degree',3,'i',1,'j',1), ...
            'CR3_x2y1',struct('degree',3,'i',2,'j',1), ...
            'CR3_y3',struct('degree',3,'i',0,'j',3))
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

        function testH1Error(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,fe] = exEllipticCRAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,2*fe.degree,'exact',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testL2Error(testCase,monom)
            p = problemScalarEllipticSquareMonomExact(monom.i,monom.j);
            [x,fe] = exEllipticCRAdaptivePatchTest.solve(testCase.mesh,p,monom.degree);

            err = ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe] = solve(mesh,p,degree)
            fe = CRk(degree,mesh);
            A = BilinearForm('alpha(x)grad(u)*grad(v)dx + beta(x)*grad(u)vdx + gamma(x)u*v dx');
            A = A.assemble(fe,'alpha',@p.alpha,'beta',@p.beta,'gamma',@p.gamma);
            b = LinearForm('f(x) v dx + g(x) v dsn').assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
            x = fe.interpolateBoundary(@p.u4Db);
            b = b - A*x;
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
        end

    end

end
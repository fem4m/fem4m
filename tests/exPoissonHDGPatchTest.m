classdef exPoissonHDGPatchTest < PatchTestCase
    %exPoissonHDGPatchTest - Patch test of the example exPoissonHDG
    %   The hybridizable dG method with flux space DGDkVector(degree),
    %   potential space DGDk(degreeU), hybrid space
    %   DGDkEdge(max(degree,degreeU)), stabilization tau = 1 and
    %   homogeneous Dirichlet data reproduces the biquadratic solution of
    %   problemScalarEllipticSquareExact for degreeU >= 4 and
    %   degree >= 3, both with static condensation ('sc') and with the
    %   monolithic system.
    %
    %   See also exPoissonHDG, problemScalarEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degrees = struct('DGD3_DGD4',struct('degree',3,'degreeU',4), ...
            'DGD4_DGD4',struct('degree',4,'degreeU',4))
        solver = struct('sc','sc','monolithic','monolithic')
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

        function testL2Error(testCase,degrees,solver)
            p = problemScalarEllipticSquareExact();
            [x,fe,dg] = exPoissonHDGPatchTest.solve(testCase.mesh,p,degrees,solver);

            err = ErrorForm('(exact(x)-u)^2dx').assemble(x(fe.nrDoF+(1:dg.nrDoF)),dg, ...
                2*dg.degree+2,'exact',@p.exactU);

            testCase.verifyLessThan(err,testCase.tol);
        end

        function testSigmaError(testCase,degrees,solver)
            p = problemScalarEllipticSquareExact();
            [x,fe] = exPoissonHDGPatchTest.solve(testCase.mesh,p,degrees,solver);

            err = ErrorForm('(exactGrad(x)-u)^2dx').assemble(-x(1:fe.nrDoF),fe, ...
                2*fe.degree+2,'exactGrad',@p.exactGradU);

            testCase.verifyLessThan(err,testCase.tol);
        end

    end

    methods (Static)

        function [x,fe,dg] = solve(mesh,p,degrees,solver)
            tau = 1;
            fe = DGDkVector(degrees.degree,mesh);
            dg = DGDk(degrees.degreeU,mesh);
            dgE = DGDkEdge(max(degrees.degreeU,degrees.degree),mesh);
            A = BilinearForm('sigma*tau dx',{'sigma','tau'}).assemble(fe);
            B = BilinearForm('div(sigma) v dx',{'sigma','v'}).assemble(fe,dg);
            C = BilinearForm('[sigma]*n {mu} dsi',{'sigma','mu'}).assemble(fe,dgE);
            C = C(dgE.freeNodes,:);
            D1 = BilinearForm('t(x)*{u}{v} dsi',{'u','v'}).assemble(dg,dg,'t',tau);
            D2 = BilinearForm('t(x)*[u][v] dsi',{'u','v'}).assemble(dg,dg,'t',tau);
            DBD = BilinearForm('t(x)*{u}{v} dsd',{'u','v'}).assemble(dg,dg,'t',tau);
            E = BilinearForm('t(x)*{u}{lambda} dsi',{'u','lambda'}).assemble(dg,dgE,'t',tau);
            E = E(dgE.freeNodes,:);
            F = BilinearForm('t(x)*{lambda}{mu} dsi',{'lambda','mu'}).assemble(dgE,dgE,'t',tau)';
            F = F(dgE.freeNodes,dgE.freeNodes);
            b = LinearForm('f(x)vdx').assemble(dg,'f',@p.f);
            if strcmpi(solver,'sc')
                % static condensation as in the example, A11 is block
                % diagonal with one block per triangle
                A11 = [A, -B'; B, 2*D1+1/2.*D2+DBD];
                A12 = [C'; -2*E'];
                A21 = [C, 2*E];
                A22 = -2*F;
                r = [zeros(fe.nrDoF,1); b];
                lambda = (A21*(A11\A12) - A22)\(A21*(A11\r));
                x = A11\(r - A12*lambda);
            else
                x = [A, -B', C';
                    B, 2*D1+1/2.*D2+DBD, -2*E';
                    C, 2*E, -2*F] \ [zeros(fe.nrDoF,1); b; zeros(numel(dgE.freeNodes),1)];
            end
        end

    end

end
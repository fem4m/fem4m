classdef ConvergenceTest < UnitTestCase
    %ConvergenceTest - Convergence rates for smooth exact solutions
    %   The Poisson problem -Delta u = f with homogeneous Dirichlet data is
    %   solved on a sequence of uniformly refined meshes for the smooth
    %   exact solutions u = sin(pi x)sin(pi y) on the unit square
    %   (problemScalarEllipticSquareOscExact) and
    %   u = sin(pi x)sin(pi y)sin(pi z) on the unit cube
    %   (problemScalarEllipticCubeExact2), and the Stokes problem for the
    %   trigonometric flow of problemStokesExactSquareSinus. The empirical
    %   order
    %       log(err(end-1)/err(end)) / log(h(end-1)/h(end))
    %   of the last refinement, with the maximal edge length h, has to be
    %   within orderTol of the theoretical order:
    %       Pk, Pk3d      |u-u_h|_H1 = O(h^k),  ||u-u_h||_L2 = O(h^(k+1))
    %       CRk           |u-u_h|_H1(T) = O(h^k), ||u-u_h||_L2 = O(h^(k+1))
    %       RTk x DGDk(k-1) mixed method
    %                     ||sigma-sigma_h||_L2 = O(h^k), ||u-u_h||_L2 = O(h^k)
    %       DGDk symmetric interior penalty (SIPDG) method with penalty
    %       2k(k+1)/h
    %                     |u-u_h|_H1(T) = O(h^k), ||u-u_h||_L2 = O(h^(k+1))
    %       Taylor-Hood PkVector x Pk(k-1) with mean-zero pressure
    %                     |u-u_h|_H1 = O(h^k),  ||u-u_h||_L2 = O(h^(k+1)),
    %                     ||p-p_h||_L2 = O(h^k)
    %   The double eigenvalue 5 pi^2 (mode 2) of the twisted diffusion
    %   eigenvalue problem of problemScalarEllipticEigenvalueSquareTwist
    %   is approximated with Pk, and the errors of the eigenfunction are
    %   the distances to the exact eigenspace:
    %       Pk eigenvalue |u*-u_h|_H1 = O(h^k), ||u*-u_h||_L2 = O(h^(k+1)),
    %                     |lambda-lambda_h| = O(h^(2k)), u* = P_E u_h
    %   The observed orders are logged at verbosity level Detailed.
    %
    %   ConvergenceTest properties:
    %       orderTol  - Tolerance of the empirical orders
    %
    %   See also problemScalarEllipticSquareOscExact,
    %   problemScalarEllipticCubeExact2, problemStokesExactSquareSinus,
    %   problemScalarEllipticEigenvalueSquareTwist,
    %   Pk, Pk3d, CRk, RTk, DGDk, PkVector, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (Constant)
        orderTol = 0.2
    end

    properties (TestParameter)
        method = struct( ...
            'P1',struct('name','Pk','degree',1,'levels',5,'order',[1 2]), ...
            'P2',struct('name','Pk','degree',2,'levels',4,'order',[2 3]), ...
            'P3',struct('name','Pk','degree',3,'levels',3,'order',[3 4]), ...
            'CR1',struct('name','CRk','degree',1,'levels',4,'order',[1 2]), ...
            'CR3',struct('name','CRk','degree',3,'levels',3,'order',[3 4]), ...
            'RT1',struct('name','RTk','degree',1,'levels',4,'order',[1 1]), ...
            'RT2',struct('name','RTk','degree',2,'levels',3,'order',[2 2]), ...
            'RT3',struct('name','RTk','degree',3,'levels',3,'order',[3 3]), ...
            'SIPDG1',struct('name','DGDk','degree',1,'levels',5,'order',[1 2]), ...
            'SIPDG2',struct('name','DGDk','degree',2,'levels',4,'order',[2 3]), ...
            'SIPDG3',struct('name','DGDk','degree',3,'levels',3,'order',[3 4]), ...
            'TH2',struct('name','PkVector','degree',2,'levels',4,'order',[2 3 2]), ...
            'TH3',struct('name','PkVector','degree',3,'levels',4,'order',[3 4 3]), ...
            'P1_3d',struct('name','Pk3d','degree',1,'levels',5,'order',[1 2]), ...
            'P2_3d',struct('name','Pk3d','degree',2,'levels',4,'order',[2 3]), ...
            'P1_eig',struct('name','PkEigenvalue','degree',1,'levels',5,'order',[1 2 2]), ...
            'P2_eig',struct('name','PkEigenvalue','degree',2,'levels',4,'order',[2 3 4]))
    end

    methods (Test)

        function testOrder(testCase,method)
            [order,err,h] = ConvergenceTest.orders(method);

            names = ConvergenceTest.errorNames(method.name);
            msg = sprintf('%s(%d): h = %s',method.name,method.degree,mat2str(h,3));
            for j = 1 : numel(method.order)
                msg = [msg,sprintf(', %s = %s (order %.3f)',names{j}, ...
                    mat2str(err(:,j)',3),order(j))]; %#ok<AGROW>
            end
            testCase.log(matlab.unittest.Verbosity.Detailed,msg);
            for j = 1 : numel(method.order)
                testCase.verifyEqual(order(j),method.order(j), ...
                    'AbsTol',ConvergenceTest.orderTol, ...
                    sprintf('empirical order of the %s',names{j}));
            end
        end

    end

    methods (Static)

        function [order,err,h] = orders(method)
            %orders - Empirical orders of the last uniform refinement
            %   [ORDER,ERR,H] = ConvergenceTest.orders(METHOD) solves on
            %   METHOD.levels uniformly refined meshes and returns the
            %   empirical orders ORDER of the last refinement, the errors
            %   ERR (one row per level, one column per error) and the mesh
            %   sizes H.
            switch method.name
                case 'Pk3d'
                    p = problemScalarEllipticCubeExact2();
                    mesh = Mesh3d(p.geom);
                case 'PkVector'
                    p = problemStokesExactSquareSinus();
                    mesh = Mesh(p.geom);
                case 'PkEigenvalue'
                    % the twist needs a few elements to be resolved
                    p = problemScalarEllipticEigenvalueSquareTwist();
                    mesh = Mesh(p.geom).refine('uniform',2);
                otherwise
                    p = problemScalarEllipticSquareOscExact();
                    mesh = Mesh(p.geom);
            end
            err = zeros(method.levels,numel(method.order));
            h = zeros(1,method.levels);
            for level = 1 : method.levels
                mesh = mesh.refine('uniform');
                h(level) = max(mesh.length4ed);
                switch method.name
                    case 'RTk'
                        err(level,:) = ConvergenceTest.solveMixed(p,mesh,method.degree);
                    case 'DGDk'
                        err(level,:) = ConvergenceTest.solveSIPDG(p,mesh,method.degree);
                    case 'PkVector'
                        err(level,:) = ConvergenceTest.solveTaylorHood(p,mesh,method.degree);
                    case 'PkEigenvalue'
                        err(level,:) = ConvergenceTest.solveEigenvalue(p,mesh,method.degree,2);
                    otherwise
                        err(level,:) = ConvergenceTest.solve(p,mesh,method);
                end
            end
            order = log(err(end-1,:)./err(end,:))/log(h(end-1)/h(end));
        end

        function err = solve(p,mesh,method)
            %solve - H1 and L2 errors of the primal Poisson discretization
            fe = feval(method.name,method.degree,mesh);
            degree = 2*method.degree + 2;
            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);
            b = LinearForm('f(x) v dx').assemble(fe,degree,'f',@p.f);
            x = zeros(fe.nrDoF,1);
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);
            err = ConvergenceTest.errors(p,x,fe,degree);
        end

        function err = solveSIPDG(p,mesh,degree)
            %solveSIPDG - H1 and L2 errors of the SIPDG method with DGDk
            % the homogeneous Dirichlet data is imposed weakly by the
            % boundary terms of 'ds'
            fe = DGDk(degree,mesh);
            quadDegree = 2*degree + 2;
            a = BilinearForm(['grad(u)*grad(v)dx', ...
                '- [u]*n.{grad(v)}ds - {grad(u)}.[v]*nds', ...
                '+ h^-1*alpha(x)*[u][v] ds']);
            A = a.assemble(fe,'alpha',2*degree*(degree+1));
            b = LinearForm('f(x) v dx').assemble(fe,quadDegree,'f',@p.f);
            x = A\b;
            err = ConvergenceTest.errors(p,x,fe,quadDegree);
        end

        function err = solveMixed(p,mesh,degree)
            %solveMixed - Flux and L2 errors of the mixed RTk x DGDk(k-1) method
            % the Dirichlet data is natural, all degrees of freedom are free
            rt = RTk(degree,mesh);
            dg = DGDk(degree-1,mesh);
            quadDegree = 2*degree + 2;
            A = BilinearForm('sigma*tau dx',{'sigma','tau'}).assemble(rt);
            B = BilinearForm('div(sigma)*v dx',{'sigma','v'}).assemble(rt,dg);
            b = LinearForm('f(x)*v dx').assemble(dg,quadDegree,'f',@(z)(-p.f(z)));
            uD = LinearForm('bd(x)*tau.n dsd',{'sigma','tau'}).assemble(rt, ...
                quadDegree,'bd',@p.u4Db);
            x = [A,B';B,sparse(dg.nrDoF,dg.nrDoF)]\[uD;b];
            err = [ErrorForm('(exact(x)-sigma)^2 dx','sigma').assemble( ...
                x(1:rt.nrDoF),rt,quadDegree,'exact',@p.exactSigma), ...
                ErrorForm('(exact(x)-u)^2 dx').assemble( ...
                x(rt.nrDoF+1:end),dg,quadDegree,'exact',@p.exactU)];
        end

        function err = solveTaylorHood(p,mesh,degree)
            %solveTaylorHood - Errors of the Taylor-Hood Stokes discretization
            %   ERR = [|u-u_h|_H1, ||u-u_h||_L2, ||p-p_h||_L2] for the
            %   velocity in PkVector(k) with Dirichlet data interpolated in
            %   the Lagrange points and the pressure in Pk(k-1) with mean
            %   zero imposed by a Lagrange multiplier.
            feV = PkVector(degree,mesh);
            feP = Pk(degree-1,mesh);
            quadDegree = 2*degree + 2;
            A = BilinearForm('grad(u)*grad(v) dx').assemble(feV);
            B = BilinearForm('-div(u)*v dx').assemble(feV,feP);
            C = LinearForm('v dx').assemble(feP);
            b = LinearForm('f(x)*v dx').assemble(feV,quadDegree,'f',@p.f);
            x = zeros(feV.nrDoF+feP.nrDoF,1);
            fixed = feV.fe{1}.fixedNodes;
            x(feV.fixedNodes) = reshape(p.u4Db(feV.fe{1}.c4dof(fixed,:)),[],1);
            b = b - A*x(1:feV.nrDoF);
            free = feV.freeNodes(:);
            S = [A(free,free),B(:,free)',sparse(numel(free),1);
                B(:,free),sparse(feP.nrDoF,feP.nrDoF),C;
                sparse(1,numel(free)),C',0];
            y = S\[b(free);-B*x(1:feV.nrDoF);0];
            x([free;feV.nrDoF+(1:feP.nrDoF)']) = y(1:end-1);
            u = x(1:feV.nrDoF);
            err = [ErrorForm('(exact(x)-grad(u))^2 dx').assemble(u,feV,quadDegree, ...
                'exact',@p.exactGradU), ...
                ErrorForm('(exact(x)-u)^2 dx').assemble(u,feV,quadDegree, ...
                'exact',@p.exactU), ...
                ErrorForm('(exact(x)-u)^2 dx').assemble(x(feV.nrDoF+1:end),feP, ...
                quadDegree,'exact',@p.exactP)];
        end

        function err = solveEigenvalue(p,mesh,degree,m)
            %solveEigenvalue - Errors of the Pk diffusion eigenvalue problem
            %   ERR = [|u*-u_h|_H1, ||u*-u_h||_L2, |lambda-lambda_h|] for
            %   mode M with ||u_h||_L2 = 1 and the L2-projection u* of u_h
            %   onto the exact eigenspace, which fixes the sign and the
            %   limit of u_h within a multiple eigenspace.
            fe = Pk(degree,mesh);
            quadDegree = 2*degree + 2;
            free = fe.freeNodes;
            A = BilinearForm('alpha(x)grad(u)*grad(v)dx').assemble(fe,'alpha',@p.alpha);
            B = BilinearForm('u*v dx').assemble(fe);
            [V,D] = eigs(A(free,free),B(free,free),max(p.cluster(m)),'smallestabs');
            [D,ind] = sort(diag(D));
            x = zeros(fe.nrDoF,1);
            x(free) = V(:,ind(m));
            x = x/sqrt(x'*B*x);
            pairs4l = p.eigenspace(m);
            uStar = @(z)0;
            gradUStar = @(z)0;
            for l = 1 : size(pairs4l,1)
                j = pairs4l(l,1);
                k = pairs4l(l,2);
                a = LinearForm('f(x) v dx').assemble(fe,quadDegree,'f',@(z)p.exactU(z,j,k))'*x;
                uStar = @(z)uStar(z) + a*p.exactU(z,j,k);
                gradUStar = @(z)gradUStar(z) + a*p.exactGradU(z,j,k);
            end
            err = [ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,quadDegree, ...
                'exact',gradUStar), ...
                ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,quadDegree, ...
                'exact',uStar), ...
                abs(D(m) - p.eigenvalue(m))];
        end

        function err = errors(p,x,fe,degree)
            %errors - (Broken) H1 seminorm and L2 norm of the error
            err = [ErrorForm('( exact(x) - grad(u) )^2 dx').assemble(x,fe,degree, ...
                'exact',@p.exactGradU), ...
                ErrorForm('( exact(x) - u )^2 dx').assemble(x,fe,degree, ...
                'exact',@p.exactU)];
        end

        function names = errorNames(name)
            %errorNames - Names of the errors of a method
            switch name
                case 'RTk'
                    names = {'L2 error of sigma','L2 error of u'};
                case 'PkVector'
                    names = {'H1 error of u','L2 error of u','L2 error of p'};
                case 'PkEigenvalue'
                    names = {'H1 error of u','L2 error of u','eigenvalue error'};
                otherwise
                    names = {'H1 error','L2 error'};
            end
        end

    end

end
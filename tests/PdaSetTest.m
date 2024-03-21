classdef PdaSetTest < UnitTestCase
    %PdaSetTest - Unit tests of the primal-dual active set method
    %   pdaSet solves the discrete obstacle problem of the 1d Laplacian
    %   with a lower obstacle, an upper obstacle, or both. The solution
    %   satisfies the complementarity conditions: the multipliers
    %   LAMBDA = A*X - B vanish where the bounds are not attained, are
    %   nonnegative on the lower and nonpositive on the upper contact set,
    %   and the bounds hold. Only the entries ISVERTEX are constrained.
    %   The inner solves with a preconditioned GMRES give the same
    %   solution as the direct solves, and the diagnostic output reports
    %   the iterations and a stop without convergence.
    %
    %   See also pdaSet, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (Constant)
        n = 39
    end

    properties (TestParameter)
        obstacle = struct( ...
            'lower',struct('f',-10,'low',-0.2,'up',[]), ...
            'upper',struct('f',10,'low',[],'up',0.2), ...
            'both',struct('f','sin','low',-0.05,'up',0.05))
    end

    methods (Test)

        function testObstacle(testCase,obstacle)
            [A,b,low,up] = PdaSetTest.problem(obstacle);

            [x,lambda,k] = pdaSet(A,b,low,up);

            PdaSetTest.verifyComplementarity(testCase,A,b,low,up,x,lambda,true(size(b)));
            testCase.verifyGreaterThan(k,0,'active set iterations');
        end

        function testVertexConstraints(testCase)
            [A,b,low,up] = PdaSetTest.problem(testCase.obstacle.both);
            isVertex = mod(1:numel(b),2)' == 1;

            [x,lambda] = pdaSet(A,b,low,up,[],[],isVertex);

            PdaSetTest.verifyComplementarity(testCase,A,b,low,up,x,lambda,isVertex);
            testCase.verifyEqual(lambda(~isVertex),zeros(nnz(~isVertex),1),'unconstrained entries');
        end

        function testPreconditionedSolver(testCase)
            [A,b,low,up] = PdaSetTest.problem(testCase.obstacle.both);
            expected = pdaSet(A,b,low,up);
            P = @(r)A\r;

            [x,lambda,~,iter] = pdaSet(A,b,low,up,[],[],[],0,P,1e-12,1e-14);

            testCase.verifyEqual(x,expected,'AbsTol',1e-8,'same solution');
            testCase.verifyGreaterThan(iter,0,'inner GMRES iterations');
            PdaSetTest.verifyComplementarity(testCase,A,b,low,up,x,lambda,true(size(b)));
        end

        function testInitialIterate(testCase)
            [A,b,low,up] = PdaSetTest.problem(testCase.obstacle.lower);
            expected = pdaSet(A,b,low,up);

            x = pdaSet(A,b,low,up,[],low);

            testCase.verifyEqual(x,expected,'AbsTol',testCase.tol);
        end

        function testNoConvergence(testCase)
            [A,b,low,up] = PdaSetTest.problem(testCase.obstacle.both); %#ok<ASGLU> used in evalc

            output = evalc('pdaSet(A,b,low,up,1);');

            testCase.verifySubstring(output,'Stopped without converging');
        end

        function testDiagnosticOutput(testCase)
            [A,b,low,up] = PdaSetTest.problem(testCase.obstacle.both); %#ok<ASGLU> used in evalc
            P = @(r)A\r; %#ok<NASGU> used in evalc

            counter = evalc('pdaSet(A,b,low,up,[],[],[],1);');
            details = evalc('pdaSet(A,b,low,up,[],[],[],3);');
            preconditioned = evalc('pdaSet(A,b,low,up,[],[],[],3,P);');

            testCase.verifyNotEmpty(counter,'iteration counter');
            testCase.verifySubstring(details,'global res:');
            testCase.verifySubstring(details,'Iterations =');
            testCase.verifySubstring(preconditioned,'Iterations =');
        end

    end

    methods (Static)

        function [A,b,low,up] = problem(obstacle)
            %problem - Finite differences of -u'' = f on (0,1), u(0)=u(1)=0
            n = PdaSetTest.n;
            h = 1/(n+1);
            x = h*(1:n)';
            e = ones(n,1);
            A = spdiags([-e 2*e -e],-1:1,n,n)/h;
            if ischar(obstacle.f)
                f = 20*sin(2*pi*x);
            else
                f = obstacle.f*e;
            end
            b = h*f;
            low = [];
            up = [];
            if ~isempty(obstacle.low)
                low = obstacle.low*e;
            end
            if ~isempty(obstacle.up)
                up = obstacle.up*e;
            end
        end

        function verifyComplementarity(testCase,A,b,low,up,x,lambda,isVertex)
            %verifyComplementarity - Complementarity of the solution
            %   Empty bounds are infinite.
            tol = 1e-8;
            if isempty(low)
                low = -inf(size(b));
            end
            if isempty(up)
                up = inf(size(b));
            end
            testCase.verifyEqual(lambda,A*x - b,'AbsTol',tol,'multipliers');
            lower = lambda > tol;
            upper = lambda < -tol;
            testCase.assertTrue(any(lower | upper),'the obstacle is attained');
            testCase.verifyEqual(x(lower),low(lower),'AbsTol',tol,'lower contact');
            testCase.verifyEqual(x(upper),up(upper),'AbsTol',tol,'upper contact');
            testCase.verifyGreaterThanOrEqual(x(isVertex) - low(isVertex),-tol,'lower bound');
            testCase.verifyLessThanOrEqual(x(isVertex) - up(isVertex),tol,'upper bound');
        end

    end

end
classdef FormEvalTest < UnitTestCase
    %FormEvalTest - Unit tests of the coefficient and function evaluation
    %   FormEval.evalCoefficient turns a coefficient into array data that
    %   acts on values of a FormType: a scalar for scalar and tensor, a
    %   vector mapping a vector type to scalars, or an n x n matrix for the
    %   n components of a vector or matrix type on every point. The tests
    %   cover every kind of coefficient, empty, scalar, constant vector
    %   or matrix, one value per element, and function handles of the
    %   points or of the points, reference points and parts, for every
    %   shape: scalar coefficients give the scalar times the identity of
    %   the shape, full coefficients are passed on. FormEval.evalFunction
    %   is tested for every signature of the source function.
    %
    %   See also FormEval, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (Constant)
        pts = [0.1 0.2;0.3 0.5;0.7 0.4]
        pts_ref = [0.2 0.3]
        parts = [2 4 5]
    end

    properties (TestParameter)
        types = struct('scalar',{{FormType.scalar}}, ...
            'vector',{{FormType.vector2d,FormType.scalar}}, ...
            'vector2d',{{FormType.vector2d}},'vector3d',{{FormType.vector3d}}, ...
            'symmatrix2d',{{FormType.symmatrix2d}},'matrix2d',{{FormType.matrix2d}}, ...
            'symmatrix3d',{{FormType.symmatrix3d}},'matrix3d',{{FormType.matrix3d}}, ...
            'tensor',{{FormType.tensor}})
        matrixTypes = struct('vector',{{FormType.vector2d,FormType.scalar}}, ...
            'vector2d',{{FormType.vector2d}},'vector3d',{{FormType.vector3d}}, ...
            'symmatrix2d',{{FormType.symmatrix2d}},'matrix2d',{{FormType.matrix2d}}, ...
            'symmatrix3d',{{FormType.symmatrix3d}},'matrix3d',{{FormType.matrix3d}})
    end

    methods (Test)

        function testEmptyCoefficient(testCase,types)
            val = FormEval.evalCoefficient(testCase.pts,testCase.pts_ref,testCase.parts,[],types{:});

            testCase.verifyEqual(val,FormEvalTest.identity(types));
        end

        function testScalarCoefficient(testCase,types)
            val = FormEval.evalCoefficient(testCase.pts,testCase.pts_ref,testCase.parts,2.5,types{:});

            testCase.verifyEqual(val,2.5*FormEvalTest.identity(types));
        end

        function testConstantCoefficient(testCase,matrixTypes)
            n = FormEvalTest.size4types(matrixTypes);
            if FormEvalTest.isVectorCoefficient(matrixTypes)
                alpha = [3 4];
                expected = alpha;
            else
                alpha = reshape(1:n^2,n,n);
                expected = reshape(alpha,[1 n n]);
            end

            val = FormEval.evalCoefficient(testCase.pts,testCase.pts_ref,testCase.parts,alpha,matrixTypes{:});

            testCase.verifyEqual(val,expected);
        end

        function testElementwiseCoefficient(testCase,types)
            alpha = (11:16)';

            val = FormEval.evalCoefficient(testCase.pts,testCase.pts_ref,testCase.parts,alpha,types{:});

            testCase.verifyEqual(val,alpha(testCase.parts).*FormEvalTest.identity(types));
        end

        function testScalarFunction(testCase,types)
            f = @(x)x(:,1) + 2*x(:,2);
            g = @(x,x_ref,parts)x(:,1) + parts(:);
            expected = FormEvalTest.identity(types);
            if numel(expected) == 1 || FormEvalTest.isVectorCoefficient(types)
                expected = 1;
            end

            val = FormEval.evalCoefficient(testCase.pts,testCase.pts_ref,testCase.parts,f,types{:});
            testCase.verifyEqual(val,f(testCase.pts).*expected,'AbsTol',testCase.tol,'f(x)');

            val = FormEval.evalCoefficient(testCase.pts,testCase.pts_ref,testCase.parts,g,types{:});
            testCase.verifyEqual(val,g(testCase.pts,[],testCase.parts).*expected, ...
                'AbsTol',testCase.tol,'f(x,x_ref,parts)');
        end

        function testFullFunction(testCase,matrixTypes)
            n = FormEvalTest.size4types(matrixTypes);
            if FormEvalTest.isVectorCoefficient(matrixTypes)
                full = @(x)x(:,1)*[1 2] + x(:,2);
            else
                full = @(x)x(:,1).*reshape(1:n^2,[1 n n]) + x(:,2);
            end
            g = @(x,x_ref,parts)full(x) + parts(:);

            val = FormEval.evalCoefficient(testCase.pts,testCase.pts_ref,testCase.parts,full,matrixTypes{:});
            testCase.verifyEqual(val,full(testCase.pts),'AbsTol',testCase.tol,'f(x)');

            val = FormEval.evalCoefficient(testCase.pts,testCase.pts_ref,testCase.parts,g,matrixTypes{:});
            testCase.verifyEqual(val,g(testCase.pts,[],testCase.parts), ...
                'AbsTol',testCase.tol,'f(x,x_ref,parts)');
        end

        function testInvalidCoefficient(testCase)
            testCase.verifyError(@()FormEval.evalCoefficient(testCase.pts, ...
                testCase.pts_ref,testCase.parts,[1 2],FormType.vector2d),'fem4m:FormEval:invalidCoefficient');
            testCase.verifyError(@()FormEval.evalCoefficient(testCase.pts, ...
                testCase.pts_ref,testCase.parts,{1},FormType.scalar),'fem4m:FormEval:invalidCoefficient');
            testCase.verifyError(@()FormEval.evalCoefficient(testCase.pts, ...
                testCase.pts_ref,testCase.parts,1,3),'fem4m:FormEval:invalidType');
            testCase.verifyError(@()FormEval.evalCoefficient(testCase.pts, ...
                testCase.pts_ref,testCase.parts,1,FormType.matrix2d,FormType.scalar),'fem4m:FormEval:invalidType');
        end

        function testEvalFunction(testCase)
            x = testCase.pts;
            parts = testCase.parts;
            normals = [1 0;0 1;-1 0];

            testCase.verifyEqual(FormEval.evalFunction(x,[],parts,[]),zeros(3,1),'empty');
            testCase.verifyEqual(FormEval.evalFunction(x,[],parts,2),2*ones(3,1),'constant');
            testCase.verifyEqual(FormEval.evalFunction(x,[],parts,@(x)x(:,1)), ...
                x(:,1),'f(x)');
            testCase.verifyEqual(FormEval.evalFunction(x,[],parts,@(x,n)sum(x.*n,2),normals), ...
                sum(x.*normals,2),'f(x,n)');
            testCase.verifyEqual(FormEval.evalFunction(x,[],parts,@(x,x_ref,parts)parts(:)), ...
                parts(:),'f(x,x_ref,parts)');
            testCase.verifyEqual(FormEval.evalFunction(x,[],parts, ...
                @(x,x_ref,parts,n)parts(:).*n(:,1),normals),parts(:).*normals(:,1),'f(x,x_ref,parts,n)');
            testCase.verifyEqual(FormEval.evalFunction(x,[],parts,@(varargin)varargin{1}(:,2)), ...
                x(:,2),'f(varargin)');
        end

    end

    methods (Static)

        function n = size4types(types)
            %size4types - Number of columns of the coefficient for TYPES
            %   TYPES is the cell {TYPE} or {TYPE,TARGET} of FormType
            %   arguments; tensor coefficients are scalars.
            n = types{1}.nrComponents;
            if isnan(n)
                n = 1;
            end
        end

        function tf = isVectorCoefficient(types)
            %isVectorCoefficient - True for vector to scalar coefficients
            tf = numel(types) == 2 && types{2} == 'scalar' && types{1}.isVector();
        end

        function val = identity(types)
            %identity - Identity of the coefficient for TYPES
            %   Vector coefficients are ones, the others the 1 x n x n
            %   identity.
            n = FormEvalTest.size4types(types);
            if FormEvalTest.isVectorCoefficient(types)
                val = ones(1,n);
            else
                val = reshape(eye(n),[1 n n]);
            end
        end

    end

end
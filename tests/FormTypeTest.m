classdef FormTypeTest < UnitTestCase
    %FormTypeTest - Unit tests of the value types of the form classes
    %   FormType.of classifies values by their number of components and
    %   the space dimension; the three components of a 2d symmetric
    %   matrix and of a 3d vector are told apart by the dimension, and
    %   component counts or arrays that match no type give tensor. The
    %   index data of every member is checked against explicit matrices:
    %   full expands the stored components to the full matrix, diagonal
    %   picks the trace entries and toMatrix reshapes into dim x dim
    %   blocks along any array dimension. symPart gives the stored
    %   components of the symmetric part of a matrix, along any array
    %   dimension and with trailing dimensions, and keeps symmetric
    %   matrices.
    %
    %   See also FormType, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        matrixType = struct('matrix2d',FormType.matrix2d,'matrix3d',FormType.matrix3d, ...
            'symmatrix2d',FormType.symmatrix2d,'symmatrix3d',FormType.symmatrix3d)
    end

    methods (Test)

        function testOf(testCase)
            testCase.verifyEqual(FormType.of(1,2),FormType.scalar);
            testCase.verifyEqual(FormType.of(1,3),FormType.scalar);
            testCase.verifyEqual(FormType.of(2,2),FormType.vector2d);
            testCase.verifyEqual(FormType.of(3,2),FormType.symmatrix2d);
            testCase.verifyEqual(FormType.of(4,2),FormType.matrix2d);
            testCase.verifyEqual(FormType.of(3,3),FormType.vector3d);
            testCase.verifyEqual(FormType.of(6,3),FormType.symmatrix3d);
            testCase.verifyEqual(FormType.of(9,3),FormType.matrix3d);
        end

        function testOfTensor(testCase)
            testCase.verifyEqual(FormType.of(5,2),FormType.tensor);
            testCase.verifyEqual(FormType.of(4,3),FormType.tensor);
            testCase.verifyEqual(FormType.of([4 2],2),FormType.tensor);
            testCase.verifyEqual(FormType.of([3 1],3),FormType.vector3d);
        end

        function testKinds(testCase)
            testCase.verifyTrue(FormType.scalar.isScalar());
            testCase.verifyTrue(FormType.vector3d.isVector());
            testCase.verifyTrue(FormType.symmatrix2d.isMatrix());
            testCase.verifyTrue(FormType.symmatrix2d.isSymmetric());
            testCase.verifyFalse(FormType.matrix2d.isSymmetric());
            testCase.verifyFalse(FormType.tensor.isMatrix());
            testCase.verifyTrue(FormType.of(3,2) == 'symmatrix2d');
        end

        function testFullAndDiagonal(testCase,matrixType)
            [A,stored] = FormTypeTest.sample(matrixType);

            testCase.verifyEqual(stored(matrixType.full),A(:)','full');
            testCase.verifyEqual(sum(stored(matrixType.diagonal)),trace(A),'diagonal');
            testCase.verifyEqual(numel(stored),matrixType.nrComponents,'nrComponents');
        end

        function testToMatrix(testCase,matrixType)
            [A,stored] = FormTypeTest.sample(matrixType);
            val = permute(stored,[1 3 2]).*[1;2];

            M = matrixType.toMatrix(val,3);

            testCase.verifyEqual(size(M),[2 1 matrixType.dim matrixType.dim]);
            testCase.verifyEqual(squeeze(M(2,1,:,:)),2*A);
            testCase.verifyEqual(matrixType.toMatrix(stored,2),reshape(A,[1 size(A)]),'dim 2');
        end

        function testToMatrixError(testCase)
            testCase.verifyError(@()FormType.vector2d.toMatrix(ones(1,2),2), ...
                'fem4m:FormType:notMatrix');
        end

        function testSymPart(testCase,matrixType)
            [A,stored] = FormTypeTest.sample(matrixType);
            S = (A + A')/2;
            if matrixType.dim == 2
                expected = [S(1,1) S(2,2) S(1,2)];
            else
                expected = [S(1,1) S(2,2) S(3,3) S(1,2) S(1,3) S(2,3)];
            end
            % rows x 1 x components x 2 sides, the second side tripled
            val = cat(4,1,3).*permute(stored,[1 3 2]).*[1;2];

            testCase.verifyEqual(matrixType.symPart(stored,2),expected,'dim 2');
            testCase.verifyEqual(matrixType.symPart(val,3), ...
                cat(4,1,3).*permute(expected,[1 3 2]).*[1;2],'two sides along dim 3');
        end

        function testSymPartError(testCase)
            testCase.verifyError(@()FormType.vector3d.symPart(ones(1,3),2), ...
                'fem4m:FormType:notMatrix');
        end

    end

    methods (Static)

        function [A,stored] = sample(type)
            %sample - Matrix of TYPE and its stored components
            n = type.dim;
            A = reshape(1:n^2,n,n);
            if type.isSymmetric()
                A = A + A';
                if n == 2
                    stored = [A(1,1) A(2,2) A(1,2)];
                else
                    stored = [A(1,1) A(2,2) A(3,3) A(1,2) A(1,3) A(2,3)];
                end
            else
                stored = A(:)';
            end
        end

    end

end
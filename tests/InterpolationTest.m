classdef InterpolationTest < UnitTestCase
    %InterpolationTest - Unit tests of the nodal interpolation
    %   The interpolation evaluates a function of the physical points,
    %   or of the points, reference points and element indices as in
    %   project, at the reference dof points c4refDof of the Lagrange
    %   elements Pk, DGPk,
    %   Pk3d and DGPk3d and of every component of the vector and matrix
    %   elements built from them. A polynomial field whose components
    %   have the degrees of the component spaces is reproduced at
    %   arbitrary reference points, also for components of different
    %   degrees and for the Mini element, whose bubble components are not
    %   interpolated. In the discontinuous spaces the interpolant of an
    %   elementwise polynomial with a different scaling on every element
    %   reproduces it on every element, and interpolating a finite
    %   element function by its getValue reproduces its coefficients.
    %
    %   See also FE, FEVector, Pk, DGPk, Pk3d, DGPk3d, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        element = struct( ...
            'P1',struct('space',@(m)Pk(1,m),'dim',2), ...
            'P2',struct('space',@(m)Pk(2,m),'dim',2), ...
            'P3',struct('space',@(m)Pk(3,m),'dim',2), ...
            'P4',struct('space',@(m)Pk(4,m),'dim',2), ...
            'P5',struct('space',@(m)Pk(5,m),'dim',2), ...
            'P1_3d',struct('space',@(m)Pk3d(1,m),'dim',3), ...
            'P2_3d',struct('space',@(m)Pk3d(2,m),'dim',3), ...
            'P3_3d',struct('space',@(m)Pk3d(3,m),'dim',3), ...
            'P4_3d',struct('space',@(m)Pk3d(4,m),'dim',3), ...
            'PVector2',struct('space',@(m)PkVector(2,m),'dim',2), ...
            'PVector2_3d',struct('space',@(m)PkVector(2,m),'dim',3), ...
            'FEVector_P2P1',struct('space',@(m)FEVector(Pk(2,m),Pk(1,m)),'dim',2), ...
            'Mini',struct('space',@(m)PkVector(1,m) + FEVector(Bk(3,m),Bk(3,m)),'dim',2))
        dgElement = struct( ...
            'DGP1',struct('space',@(m)DGPk(1,m),'dim',2), ...
            'DGP2',struct('space',@(m)DGPk(2,m),'dim',2), ...
            'DGP3',struct('space',@(m)DGPk(3,m),'dim',2), ...
            'DGP1_3d',struct('space',@(m)DGPk3d(1,m),'dim',3), ...
            'DGP2_3d',struct('space',@(m)DGPk3d(2,m),'dim',3), ...
            'DGPVector2',struct('space',@(m)DGPkVector(2,m),'dim',2), ...
            'DGPVector1_3d',struct('space',@(m)DGPkVector(1,m),'dim',3), ...
            'DGPMatrix1',struct('space',@(m)DGPkMatrix(1,m),'dim',2), ...
            'DGPMatrix2_3d',struct('space',@(m)DGPkMatrix(2,m),'dim',3), ...
            'DGPSymMatrix2',struct('space',@(m)DGPkSymMatrix(2,m),'dim',2), ...
            'DGPSymMatrix1_3d',struct('space',@(m)DGPkSymMatrix(1,m),'dim',3))
    end

    methods (Test)

        function testReproducesPolynomial(testCase,element)
            mesh = InterpolationTest.mesh(element.dim);
            fe = element.space(mesh);
            f = InterpolationTest.field(fe);

            x = fe.interpolate(f);

            InterpolationTest.verifyValues(testCase,fe,x,@(z)f(mesh.point4e(z)));
        end

        function testElementwisePolynomial(testCase,dgElement)
            mesh = InterpolationTest.mesh(dgElement.dim);
            fe = dgElement.space(mesh);
            f = InterpolationTest.field(fe);
            scale = 1 + (1:mesh.nrElems)'/mesh.nrElems;

            x = fe.interpolate(@(pts,pts_ref,parts)scale(parts).*f(pts));

            InterpolationTest.verifyValues(testCase,fe,x,@(z)scale.*f(mesh.point4e(z)));
        end

        function testFiniteElementFunction(testCase,element)
            mesh = InterpolationTest.mesh(element.dim);
            fe = element.space(mesh);
            x = fe.interpolate(InterpolationTest.field(fe));

            y = fe.interpolate(@(pts,pts_ref,parts)fe.getValue(pts,pts_ref,parts,x));

            testCase.verifyEqual(y,x,'AbsTol',testCase.tol);
        end

    end

    methods (Static)

        function verifyValues(testCase,fe,x,g)
            %verifyValues - Verify the interpolant against g at points
            %   The coefficients X of FE are contracted with the basis at
            %   interior and vertex reference points and compared with the
            %   function G of the reference coordinates.
            dim = size(fe.mesh.c4n,2);
            refPoints = [0.2 0.3 0.1;0.6 0.1 0.2;0 0 0];
            for ref = refPoints(:,1:dim)'
                basis = fe.getBasis([],ref',[]);
                val = reshape(sum(basis.*x(fe.dof4e),2),fe.mesh.nrElems,[]);
                testCase.verifyEqual(val,g(ref'),'AbsTol',testCase.tol, ...
                    sprintf('values at the reference point %s',mat2str(ref')));
            end
        end

        function mesh = mesh(dim)
            %mesh - L-shaped 2d patch mesh or 3d cube mesh
            if dim == 2
                mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            else
                mesh = UnitTestCase.cubeMesh();
            end
        end

        function f = field(fe)
            %field - Polynomial field in the space of FE
            %   F(Z) returns the components x^k - 2*x*y^(k-1) + y + 1 +
            %   c*(x - y) in 2d and x^k - 2*y^(k-1)*z + y + 1 + c*(x - y)
            %   in 3d, c = 0,1,..., of the degrees k of the scalar
            %   component spaces as columns.
            if isa(fe,'FEVector')
                degrees = cellfun(@(c)c.degree,fe.fe(:,1))';
            else
                degrees = fe.degree;
            end
            if size(fe.mesh.c4n,2) == 2
                f = @(z)z(:,1).^degrees - 2*z(:,1).*z(:,2).^(degrees-1) + z(:,2) + 1 ...
                    + (0:numel(degrees)-1).*(z(:,1) - z(:,2));
            else
                f = @(z)z(:,1).^degrees - 2*z(:,2).^(degrees-1).*z(:,3) + z(:,2) + 1 ...
                    + (0:numel(degrees)-1).*(z(:,1) - z(:,2));
            end
        end

    end

end
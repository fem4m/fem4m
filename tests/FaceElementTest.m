classdef FaceElementTest < UnitTestCase
    %FaceElementTest - Unit tests of the finite elements on edges and faces
    %   DGDkEdge is spanned by the Legendre polynomials on every edge and
    %   DGDkFace by the Dubiner polynomials on every face of a tetrahedral
    %   mesh. Both bases are orthogonal, so that the mass matrices of the
    %   averages are diagonal with the entries |E|/(2i+1) and
    %   |F|/((2i+1)(i+j+1)), the load vector of f = 1 is the measure of the
    %   sides in the first degree of freedom, and the degrees of freedom of
    %   the Dirichlet sides are fixed. The values of the second side vanish
    %   on the boundary. FEEdgeVector stacks scalar edge elements into a
    %   vector element with block diagonal mass matrix. DGPkFace is the
    %   face trace space of Pk3d; its tangential gradient is orthogonal to
    %   the face normal and its surface curl is grad x n, and the values,
    %   gradients and curls of a function contract the bases with its
    %   coefficients, with value one and vanishing derivatives for the
    %   constant one (zero on the missing second side of boundary faces).
    %
    %   See also DGDkEdge, DGDkFace, DGPkFace, FEEdge, FEEdgeVector,
    %   FEFace, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degree = struct('k0',0,'k1',1,'k2',2,'k3',3)
        traceDegree = struct('k1',1,'k2',2,'k3',3)
    end

    methods (Test)

        function testEdgeMassMatrix(testCase,degree)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = DGDkEdge(degree,mesh);

            M = BilinearForm('{u}{v} ds').assemble(fe);
            b = LinearForm('f(x){v} ds').assemble(fe,'f',@(x)ones(size(x,1),1));

            testCase.verifyEqual(fe.nrDoF,(degree+1)*mesh.nrEdges);
            expected = mesh.length4ed./(2*(0:degree)+1);
            testCase.verifyEqual(full(M),full(sparse(fe.dof4ed(:),fe.dof4ed(:),expected(:))), ...
                'AbsTol',testCase.tol,'diagonal mass matrix');
            expected = zeros(fe.nrDoF,1);
            expected(fe.dof4ed(:,1)) = mesh.length4ed;
            testCase.verifyEqual(b,expected,'AbsTol',testCase.tol,'load vector of 1');
            testCase.verifyEqual(sort(fe.fixedNodes(:)),sort(reshape(fe.dof4ed(mesh.ed4Db,:),[],1)));
        end

        function testEdgeValue(testCase,degree)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = DGDkEdge(degree,mesh);
            c = sin(1:fe.nrDoF)';
            s = 0.3;

            val = fe.getEdgeValue([],s,1:mesh.nrEdges,c);

            legendre = DGDkEdge.getEdgeRefBasis(degree,s);
            expected = sum(c(fe.dof4ed).*legendre,2);
            boundary = mesh.e4ed(:,2) == 0;
            testCase.verifySize(val,[mesh.nrEdges,2]);
            testCase.verifyEqual(val(:,1),expected,'AbsTol',testCase.tol,'first side');
            testCase.verifyEqual(val(~boundary,2),expected(~boundary),'AbsTol',testCase.tol,'second side');
            testCase.verifyEqual(val(boundary,2),zeros(nnz(boundary),1),'no second side on the boundary');
        end

        function testEdgeVector(testCase,degree)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            scalar = DGDkEdge(degree,mesh);
            fe = FEEdgeVector(scalar,DGDkEdge(degree,mesh));
            nb = scalar.nrBasisFunc;

            basis = fe.getEdgeBasis([],0.3,1:mesh.nrEdges);
            M = BilinearForm('{u}{v} ds').assemble(fe);

            testCase.verifyEqual(fe.nrDoF,2*scalar.nrDoF);
            testCase.verifyEqual(fe.nrBasisFunc,2*nb);
            testCase.verifyEqual(fe.dof4ed,[scalar.dof4ed,scalar.nrDoF + scalar.dof4ed]);
            testCase.verifyEqual(sort(fe.fixedNodes(:)),sort([scalar.fixedNodes(:);scalar.nrDoF + scalar.fixedNodes(:)]));
            testCase.verifySize(basis,[mesh.nrEdges,2*nb,2,2]);
            scalarBasis = scalar.getEdgeBasis([],0.3,1:mesh.nrEdges);
            testCase.verifyEqual(basis(:,1:nb,1,:),permute(scalarBasis,[1 2 4 3]),'first component');
            testCase.verifyEqual(basis(:,nb+1:end,2,:),permute(scalarBasis,[1 2 4 3]),'second component');
            testCase.verifyEqual(basis(:,nb+1:end,1,:),zeros(mesh.nrEdges,nb,1,2),'block structure');
            testCase.verifyEqual(basis(:,1:nb,2,:),zeros(mesh.nrEdges,nb,1,2),'block structure');
            MScalar = BilinearForm('{u}{v} ds').assemble(scalar);
            testCase.verifyEqual(full(M),full(blkdiag(MScalar,MScalar)),'AbsTol',testCase.tol, ...
                'block diagonal mass matrix');
        end

        function testFaceMassMatrix(testCase,degree)
            mesh = UnitTestCase.cubeMesh();
            fe = DGDkFace(degree,mesh);

            M = BilinearForm('{u}{v} ds').assemble(fe);
            b = LinearForm('f(x){v} ds').assemble(fe,'f',@(x)ones(size(x,1),1));

            nb = (degree+1)*(degree+2)/2;
            testCase.verifyEqual(fe.nrBasisFunc,nb);
            testCase.verifyEqual(fe.nrDoF,nb*mesh.nrFaces);
            [i,j] = FaceElementTest.dubinerIndices(degree);
            expected = mesh.area4f./((2*i+1).*(i+j+1));
            testCase.verifyEqual(full(M),full(sparse(fe.dof4f(:),fe.dof4f(:),expected(:))), ...
                'AbsTol',testCase.tol,'diagonal mass matrix');
            expected = zeros(fe.nrDoF,1);
            expected(fe.dof4f(:,1)) = mesh.area4f;
            testCase.verifyEqual(b,expected,'AbsTol',testCase.tol,'load vector of 1');
            testCase.verifyEqual(sort(fe.fixedNodes(:)),sort(reshape(fe.dof4f(mesh.f4Db,:),[],1)));
        end

        function testFaceValue(testCase)
            mesh = UnitTestCase.cubeMesh();
            fe = DGDkFace(2,mesh);
            c = sin(1:fe.nrDoF)';
            ref = [0.2 0.3];

            val = fe.getFaceValue([],ref,1:mesh.nrFaces,c);

            expected = sum(c(fe.dof4f).*DGDkFace.getFaceRefBasis(2,ref),2);
            boundary = mesh.e4f(:,2) == 0;
            testCase.verifyEqual(val(:,1),expected,'AbsTol',testCase.tol,'first side');
            testCase.verifyEqual(val(~boundary,2),expected(~boundary),'AbsTol',testCase.tol,'second side');
            testCase.verifyEqual(val(boundary,2),zeros(nnz(boundary),1),'no second side on the boundary');
        end

        function testTraceSpace(testCase,traceDegree)
            mesh = UnitTestCase.cubeMesh();

            fe = DGPkFace(traceDegree,mesh);

            ref = [0.2 0.3];
            nb = (traceDegree+1)*(traceDegree+2)/2;
            testCase.verifyEqual(fe.nrBasisFunc,nb);
            testCase.verifyEqual(fe.nrDoF,nb*mesh.nrFaces);
            basis = fe.getFaceBasis([],ref,1:mesh.nrFaces);
            testCase.verifyEqual(sum(basis(:,:,1),2),ones(mesh.nrFaces,1),'AbsTol',testCase.tol, ...
                'partition of unity');
            grad = fe.getFaceGradBasis([],ref,1:mesh.nrFaces);
            normals = mesh.normal4f;
            testCase.verifyEqual(sum(grad(:,:,:,1).*permute(normals,[1 3 2]),3), ...
                zeros(mesh.nrFaces,nb),'AbsTol',testCase.tol,'tangential gradient');
            curl = fe.getFaceCurlBasis([],ref,1:mesh.nrFaces);
            testCase.verifyEqual(curl(:,:,:,1),cross(grad(:,:,:,1), ...
                repmat(permute(normals,[1 3 2]),1,nb,1),3),'AbsTol',testCase.tol,'grad x n');
        end

        function testFaceTraceValues(testCase,traceDegree)
            mesh = UnitTestCase.cubeMesh();
            fe = DGPkFace(traceDegree,mesh);
            ref = [0.2 0.3];
            parts = mesh.nrFaces:-2:1;
            c = sin(1:fe.nrDoF)';
            one = ones(fe.nrDoF,1);
            inner = mesh.e4f(parts,2) > 0;
            c4f = c(fe.dof4f(parts,:));

            val = fe.getFaceValue([],ref,parts,c);
            basis = fe.getFaceBasis([],ref,parts);
            testCase.verifyEqual(val,permute(sum(basis.*c4f,2),[1 3 2]), ...
                'AbsTol',testCase.tol,'values');
            testCase.verifyEqual(fe.getFaceValue([],ref,parts,one),[ones(numel(parts),1),double(inner)], ...
                'AbsTol',testCase.tol,'value of one');

            grad = fe.getFaceGradValue([],ref,parts,c);
            gradBasis = fe.getFaceGradBasis([],ref,parts);
            testCase.verifyEqual(grad,permute(sum(gradBasis.*c4f,2),[1 3 4 2]), ...
                'AbsTol',testCase.tol,'tangential gradients');
            testCase.verifyEqual(fe.getFaceGradValue([],ref,parts,one),zeros(numel(parts),3,2), ...
                'AbsTol',testCase.tol,'gradient of one');

            curl = fe.getFaceCurlValue([],ref,parts,c);
            curlBasis = fe.getFaceCurlBasis([],ref,parts);
            testCase.verifyEqual(curl,permute(sum(curlBasis.*c4f,2),[1 3 4 2]), ...
                'AbsTol',testCase.tol,'surface curls');
            testCase.verifyEqual(fe.getFaceCurlValue([],ref,parts,one),zeros(numel(parts),3,2), ...
                'AbsTol',testCase.tol,'curl of one');
        end

    end

    methods (Static)

        function [i,j] = dubinerIndices(p)
            %dubinerIndices - Indices (i,j) of the Dubiner basis functions
            [j,i] = ndgrid(0:p);
            keep = i + j <= p;
            i = i(keep)';
            j = j(keep)';
        end

    end

end
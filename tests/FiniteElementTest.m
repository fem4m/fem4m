classdef FiniteElementTest < UnitTestCase
    %FiniteElementTest - Unit tests of the finite element spaces
    %   The number of degrees of freedom of the finite elements is
    %   compared with the dimension formula of the space in terms of the
    %   numbers of nodes, edges, faces and elements of the mesh; with
    %   elementwise degrees (also 0), DGDk, Bk and CurlBRotk count the
    %   dofs of the degree of every element. degree4DoF is the polynomial degree of
    %   every basis function of the hierarchical elements DGDk, Bk,
    %   CurlBRotk, BDMk, HHJk and HHJB. The divergence of every basis
    %   function of the reduced Arnold-Winther elements AW0 and AWNC0 is a
    %   rigid motion, also on the Union Jack mesh. The
    %   Lagrange basis functions form a partition of unity and are
    %   nodal (Kronecker delta at the reference dof points). The values,
    %   derivatives, Laplacians and biharmonics of a finite element
    %   function evaluated on a subset of the elements are the
    %   corresponding rows of those on all elements, and on a mesh of one
    %   element they are the basis functions contracted with the
    %   coefficients. The same holds for the divergence and rot of the
    %   matrix-valued spaces with rows on the reference element (RTk,
    %   BDMk, Nd1k, DGPk) and on the physical polygons (PDGk). The row
    %   divergence and trace of the 3d symmetric matrix spaces DGPkSymMatrix
    %   and PDGkSymMatrix are those of a linear symmetric matrix field. The
    %   discontinuous spaces, also DGPk, DGPk3d, DGHPk and Nd1DGk3d derived
    %   from continuous ones, have dof4e as their only dof table, every dof
    %   in one element, and fix nothing. The interpolations and
    %   projections are tested in InterpolationTest and ProjectionTest.
    %
    %   See also InterpolationTest, ProjectionTest, UnitTestCase, Pk,
    %   Pk3d, CRk, DGPk, DGDk, RTk, BDMk, Nd1k, Nd2k, HHJk, FEMatrix,
    %   PDGkMatrix, FESymMatrix

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        element = struct( ...
            'P1',struct('name','Pk','degree',1), ...
            'P2',struct('name','Pk','degree',2), ...
            'P4',struct('name','Pk','degree',4), ...
            'CR1',struct('name','CRk','degree',1), ...
            'CR3',struct('name','CRk','degree',3), ...
            'DGP1',struct('name','DGPk','degree',1), ...
            'DGP3',struct('name','DGPk','degree',3), ...
            'DGD2',struct('name','DGDk','degree',2), ...
            'RT1',struct('name','RTk','degree',1), ...
            'RT3',struct('name','RTk','degree',3), ...
            'BDM1',struct('name','BDMk','degree',1), ...
            'BDM2',struct('name','BDMk','degree',2), ...
            'Nd1_2',struct('name','Nd1k','degree',2), ...
            'Nd2_2',struct('name','Nd2k','degree',2), ...
            'HHJ0',struct('name','HHJk','degree',0), ...
            'HHJ2',struct('name','HHJk','degree',2))
        element3d = struct( ...
            'P1',struct('name','Pk3d','degree',1), ...
            'P2',struct('name','Pk3d','degree',2), ...
            'P4',struct('name','Pk3d','degree',4))
        lagrange = struct('P1',1,'P2',2,'P3',3,'P5',5)
        matrixElement = struct( ...
            'RT1',struct('space',@(m)FEMatrix(RTk(1,m),RTk(1,m)),'poly',false, ...
            'names',"getDivValue"), ...
            'BDM1',struct('space',@(m)FEMatrix(BDMk(1,m),BDMk(1,m)),'poly',false, ...
            'names',"getDivValue"), ...
            'Nd1_1',struct('space',@(m)FEMatrix(Nd1k(1,m),Nd1k(1,m)),'poly',false, ...
            'names',"getRotValue"), ...
            'DGP1',struct('space',@(m)DGPkMatrix(1,m),'poly',false, ...
            'names',["getDivValue" "getRotValue"]), ...
            'PDG1',struct('space',@(m)PDGkMatrix(1,m),'poly',true, ...
            'names',["getDivValue" "getRotValue"]))
        symMatrix3d = struct( ...
            'DGP1',@(m)DGPkSymMatrix(1,m), ...
            'PDG1',@(m)PDGkSymMatrix(1,m))
        elementwise = struct( ...
            'DGD1to3',struct('name','DGDk','degree',@(n)1+mod((1:n)',3),'dim',@(d)(d+1).*(d+2)/2), ...
            'DGD0to2',struct('name','DGDk','degree',@(n)mod((1:n)',3),'dim',@(d)(d+1).*(d+2)/2), ...
            'B3to5',struct('name','Bk','degree',@(n)3+mod((1:n)',3),'dim',@(d)(d-2).*(d-1)/2), ...
            'CurlBRot2to4',struct('name','CurlBRotk','degree',@(n)2+mod((1:n)',3),'dim',@(d)d.*(d+1)/2-1))
        hierarchical = struct( ...
            'DGD3',@(m)DGDk(3,m),'B4',@(m)Bk(4,m),'CurlBRot3',@(m)CurlBRotk(3,m), ...
            'BDM3',@(m)BDMk(3,m),'HHJ0',@(m)HHJk(0,m),'HHJ1',@(m)HHJk(1,m), ...
            'HHJ3',@(m)HHJk(3,m),'HHJB1',@(m)HHJB(1,m),'HHJB3',@(m)HHJB(3,m))
        reducedAW = struct( ...
            'AW0_Lshape',struct('space',@(m)AW0(3,m),'mesh',@()ProjectionTest.mesh(2)), ...
            'AW0_UnionJack',struct('space',@(m)AW0(3,m), ...
            'mesh',@()ProjectionTest.unionJackMesh().refine('uniform')), ...
            'AWNC0_Lshape',struct('space',@(m)AWNC0(2,m),'mesh',@()ProjectionTest.mesh(2)), ...
            'AWNC0_UnionJack',struct('space',@(m)AWNC0(2,m), ...
            'mesh',@()ProjectionTest.unionJackMesh().refine('uniform')))
    end

    methods (Test)

        function testNrDoF(testCase,element)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = feval(element.name,element.degree,mesh);

            [nrDoF4ed,nrDoF4e] = FiniteElementTest.dofs4entity(element.name,element.degree);
            nrDoF4n = double(any(strcmp(element.name,{'Pk'})));

            testCase.verifyEqual(fe.nrDoF, ...
                nrDoF4n*mesh.nrNodes + nrDoF4ed*mesh.nrEdges + nrDoF4e*mesh.nrElems);
            testCase.verifyEqual(fe.nrBasisFunc,3*nrDoF4n + 3*nrDoF4ed + nrDoF4e);
            testCase.verifyEqual(sort(unique(fe.dof4e(:))),(1:fe.nrDoF)');
        end

        function testValueOnSubset(testCase,element)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = feval(element.name,element.degree,mesh);
            x = sin(1:fe.nrDoF)';
            parts = mesh.nrElems:-2:1;

            for name = FiniteElementTest.valueMethods(element.name)
                val = fe.(name)([],[0.2 0.3],parts,x);

                all = fe.(name)([],[0.2 0.3],[],x);
                testCase.verifyEqual(val(:,:),all(parts,:),'AbsTol',testCase.tol,name);
            end
        end

        function testMatrixValueOnSubset(testCase,matrixElement)
            if matrixElement.poly
                mesh = MeshPolyT(Geometry(struct('name','Lshape','type','T')));
            else
                mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            end
            fe = matrixElement.space(mesh);
            x = sin(1:fe.nrDoF)';
            parts = mesh.nrElems:-2:1;
            allParts = 1:mesh.nrElems;
            if fe.isRefElemBased
                ref = [0.2 0.3];
                pts = [];
                ptsAll = [];
            else
                ref = [];
                pts = PolygonalElementTest.interiorPoints(mesh);
                ptsAll = pts;
                pts = pts(parts,:);
            end

            for name = matrixElement.names
                val = fe.(name)(pts,ref,parts,x);

                all = fe.(name)(ptsAll,ref,allParts,x);
                testCase.verifyEqual(val(:,:),all(parts,:),'AbsTol',testCase.tol,name);
            end
        end

        function testSymMatrixValues3d(testCase,symMatrix3d)
            mesh = UnitTestCase.cubeMesh();
            fe = symMatrix3d(mesh);
            ref = [0.2 0.3 0.1];
            pts = mesh.point4e(ref);
            parts = 1:mesh.nrElems;
            x = FormValueTest.coefficients(fe,'symmatrix3d');

            F = FormValueTest.evalField('symmatrix3d',pts);
            testCase.verifyEqual(fe.getDivValue(pts,ref,parts,x),F.div, ...
                'AbsTol',testCase.tol,'row divergence');
            testCase.verifyEqual(fe.getTrValue(pts,ref,parts,x),FormBasisTest.tr(F.M), ...
                'AbsTol',testCase.tol,'trace');
        end

        function testValueOnOneElement(testCase,element)
            fe = feval(element.name,element.degree,UnitTestCase.referenceTriangle());
            x = sin(1:fe.nrDoF)';
            ref = [0.2 0.3];
            u = reshape(x(fe.dof4e),1,[]);
            contract = @(basis)reshape(sum(u.*basis,2),1,[]);

            for name = FiniteElementTest.valueMethods(element.name)
                val = fe.(name)([],ref,[],x);

                switch name
                    case "getLaplaceValue"
                        d2 = contract(fe.getD2Basis([],ref,[]));
                        expected = d2(1) + d2(2);
                    case "getBiharmonicValue"
                        d4 = contract(fe.getD4Basis([],ref,[]));
                        expected = d4(1) + d4(2) + 2*d4(3);
                    otherwise
                        expected = contract(fe.(replace(name,"Value","Basis"))([],ref,[]));
                end
                testCase.verifyEqual(reshape(val,1,[]),expected,'AbsTol',testCase.tol,name);
            end
        end

        function testNrDoF3d(testCase,element3d)
            mesh = UnitTestCase.cubeMesh();
            fe = feval(element3d.name,element3d.degree,mesh);

            k = element3d.degree;
            testCase.verifyEqual(fe.nrDoF,mesh.nrNodes + (k-1)*mesh.nrEdges ...
                + (k-1)*(k-2)/2*mesh.nrFaces + (k-1)*(k-2)*(k-3)/6*mesh.nrElems);
        end

        function testNrDoFElementwise(testCase,elementwise)
            mesh = Mesh(struct('name','Square')).refine('uniform');
            degree4e = elementwise.degree(mesh.nrElems);
            fe = feval(elementwise.name,degree4e,mesh);

            testCase.verifyEqual(fe.nrDoF,sum(elementwise.dim(degree4e)));
            testCase.verifyEqual(fe.allnrDoF,mesh.nrElems*fe.nrBasisFunc);
            testCase.verifyEqual(numel(fe.freeNodes),fe.nrDoF);
        end

        function testDegree4DoF(testCase,hierarchical)
            fe = hierarchical(UnitTestCase.referenceTriangle());

            testCase.verifyEqual(fe.degree4DoF,FiniteElementTest.basisDegrees(fe));
        end

        function testRigidMotionDivergence(testCase,reducedAW)
            mesh = reducedAW.mesh();
            fe = reducedAW.space(mesh);
            % the divergence is linear: gradient from the vertex values
            p1 = Pk(1,mesh);
            gradP1 = p1.getGradBasis([],[1 1]/3,[]);
            gradDiv = zeros(mesh.nrElems,fe.nrBasisFunc,2,2);
            meanDiv = 0;
            for j = 1 : 3
                div = fe.getDivBasis([],p1.c4refDof(j,:),[]);
                gradDiv = gradDiv + div.*permute(gradP1(:,j,:),[1 2 4 3]);
                meanDiv = meanDiv + div/3;
            end
            eps = cat(3,gradDiv(:,:,1,1),gradDiv(:,:,2,2),gradDiv(:,:,1,2) + gradDiv(:,:,2,1));
            scale = max(abs(gradDiv(:)));

            testCase.verifyTrue(all(isfinite(fe.basisCoefficients(:))),'finite basis');
            testCase.verifyEqual(fe.getDivBasis([],[1 1]/3,[]),meanDiv, ...
                'AbsTol',testCase.tol*scale,'linear divergence');
            testCase.verifyEqual(eps,zeros(size(eps)),'AbsTol',testCase.tol*scale, ...
                'eps(div) = 0');
        end

        function testPartitionOfUnity(testCase,lagrange)
            mesh = PatchTestCase.patchMesh(struct('name','Square'));
            fe = Pk(lagrange,mesh);

            basis = fe.getBasis([],[0.2 0.3],[]);

            testCase.verifyEqual(sum(basis,2),ones(mesh.nrElems,1),'AbsTol',testCase.tol);
        end

        function testNodalBasis(testCase,lagrange)
            mesh = PatchTestCase.patchMesh(struct('name','Square'));
            fe = Pk(lagrange,mesh);

            values = zeros(fe.nrBasisFunc);
            for j = 1 : fe.nrBasisFunc
                basis = fe.getBasis([],fe.c4refDof(j,:),[]);
                values(:,j) = basis(1,:)';
            end

            testCase.verifyEqual(values,eye(fe.nrBasisFunc),'AbsTol',testCase.tol);
        end

        function testDofPoints(testCase,lagrange)
            mesh = PatchTestCase.patchMesh(struct('name','Square'));
            fe = Pk(lagrange,mesh);

            for j = 1 : fe.nrBasisFunc
                testCase.verifyEqual(fe.c4dof(fe.dof4e(:,j),:), ...
                    mesh.point4e(fe.c4refDof(j,:)),'AbsTol',testCase.tol);
            end
        end

        function testDiscontinuousDoFs(testCase)
            % dof4e is the only dof table of the discontinuous spaces, also
            % of those derived from a continuous sibling: every dof belongs
            % to one element and nothing is fixed
            mesh = Mesh(struct('name','Square')).refine('uniform');
            cube = UnitTestCase.cubeMesh();
            fes = {DGDk(2,mesh),DGPk(2,mesh),DGHPk(2,mesh),DGPkVector(2,mesh), ...
                DGPkMatrix(2,mesh),DGDk3d(2,cube),DGPk3d(2,cube),Nd1DGk3d(1,cube)};
            for k = 1 : numel(fes)
                fe = fes{k};
                names = properties(fe);
                names = names(startsWith(names,'dof4') & ~strcmp(names,'dof4e'));
                for j = 1 : numel(names)
                    testCase.verifyEmpty(fe.(names{j}),[class(fe),'.',names{j}]);
                end
                testCase.verifyEqual(sort(fe.dof4e(:))',1:numel(fe.dof4e),[class(fe),'.dof4e']);
                testCase.verifyEmpty(fe.fixedNodes,[class(fe),'.fixedNodes']);
            end
        end

    end

    methods (Static)

        function names = valueMethods(name)
            %valueMethods - get...Value methods supported by an element
            switch name
                case {'Pk','DGPk','DGDk'}
                    names = ["getValue" "getGradValue" "getD2Value" "getD3Value" ...
                        "getLaplaceValue" "getBiharmonicValue"];
                case 'CRk'
                    names = ["getValue" "getGradValue"];
                case {'RTk','BDMk','HHJk'}
                    names = ["getValue" "getDivValue"];
                case 'Nd1k'
                    names = ["getValue" "getRotValue" "getCurlRotValue"];
                case 'Nd2k'
                    names = ["getValue" "getRotValue"];
            end
        end

        function degree = basisDegrees(fe)
            %basisDegrees - Polynomial degree of every basis function
            %   DEGREE(J) is the smallest D for which the J-th basis
            %   function of FE on the reference triangle is a polynomial of
            %   degree D, from least squares fits in the Dubiner basis at
            %   a lattice that is unisolvent up to the degree FE.degree+3.
            n = fe.degree + 3;
            [x,y] = meshgrid((0:n)/n);
            inside = x + y <= 1 + 1e-12;
            refPoints = [x(inside),y(inside)];
            V = [];
            for k = 1 : size(refPoints,1)
                ref = refPoints(k,:);
                basis = FormBasis.elementRows(fe.getBasis(fe.mesh.point4e(ref),ref,1),1,fe);
                V = cat(1,V,reshape(basis(1,:,:),1,fe.nrBasisFunc,[]));
            end
            degree = inf(1,fe.nrBasisFunc);
            for d = 0 : n
                D = DGDk.getElementBasis(d,refPoints);
                for j = find(isinf(degree))
                    v = reshape(V(:,j,:),size(V,1),[]);
                    if norm(v - D*(D\v),inf) < 1e-9*max(1,norm(v,inf))
                        degree(j) = d;
                    end
                end
            end
        end

        function [nrDoF4ed,nrDoF4e] = dofs4entity(name,k)
            %dofs4entity - Number of edge and interior dofs of a 2d element
            switch name
                case 'Pk'
                    nrDoF4ed = k-1;
                    nrDoF4e = (k-1)*(k-2)/2;
                case 'CRk'
                    nrDoF4ed = k;
                    nrDoF4e = (k-1)*(k-2)/2;
                case {'DGPk','DGDk'}
                    nrDoF4ed = 0;
                    nrDoF4e = (k+1)*(k+2)/2;
                case {'RTk','Nd1k'}
                    nrDoF4ed = k;
                    nrDoF4e = k*(k-1);
                case {'BDMk','Nd2k'}
                    nrDoF4ed = k+1;
                    nrDoF4e = (k-1)*(k+1);
                case 'HHJk'
                    nrDoF4ed = k+1;
                    nrDoF4e = 3*(k+1)*(k+2)/2 - 3*(k+1);
            end
        end

    end

end
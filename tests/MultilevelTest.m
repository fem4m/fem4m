classdef MultilevelTest < UnitTestCase
    %MultilevelTest - Unit tests of prolongations, multigrid and PINVIT
    %   FE.prolongation maps the coefficients of a field in the coarse
    %   space to its coefficients in the space on the refined mesh. In 2d
    %   it is tested for all elements (reference based ones and those
    %   with geometry dependent bases, which use local L2 projections)
    %   under red-green-blue and newest vertex bisection refinement of
    %   some elements, two uniform refinements, graded and graded red
    %   refinement, for a degree increase, for two refine calls, for
    %   vector-valued elements, for the enrichment of HHJk by HHJB and
    %   for Argyris elements of size 1e-6 (without warnings, also of
    %   FE.project), for element-wise degrees up to 12, on strongly
    %   graded meshes and for the area-scaled basis of
    %   FESkewSymMatrix; without a
    %   refinement round it is the identity without round-off entries;
    %   the bubble spaces Bk, B2k, CurlBk, CurlBRotk, the Mini element
    %   and HHJB have prolongations of the right size. For sums of
    %   elements (FEMatrix + CurlBRotk, the Stenberg element with
    %   FEMatrix of CurlBk, the Mini element, also with elementwise
    %   bubble degrees) it is the local L2 projection onto the whole
    %   fine sum space with zero rows for the dofs that are neither free
    %   nor fixed: exact for the BDMk and RTk sums under uniform
    %   refinement, where they are nested although their bubbles are
    %   not, and otherwise for the functions of the nested first
    %   summand. DGDkVector with elementwise degrees 0 to 2, whose
    %   components are shifted by allnrDoF, has an exact prolongation.
    %   Every element is
    %   rebuilt with the same space by feval(class(fe),fe.degree,mesh),
    %   and isAffineEquivalent is false exactly for the elements with
    %   geometry dependent bases; for the affine equivalent ones the
    %   reference matrices give the same prolongation as the local L2
    %   projections on the physical elements, unless the space is
    %   Piola-mapped and not nested. In 3d it maps the coefficients of a
    %   polynomial of degree k in Pk3d and of polynomial fields in the
    %   oriented elements RTk3d, BDMk3d, Nd1k3d, Nd2k3d and HPk3d under
    %   bisection of one tetrahedron of a mesh whose other elements stay
    %   unrefined, of all tetrahedra of a graded mesh where some of them
    %   are bisected twice in one call, and under uniform refinement.
    %   The contraction number of Multigrid is the
    %   spectral radius of the error propagation I - B*A of one V-cycle
    %   B and the damping factor the spectral radius of T*A for the
    %   smoother step x - T*(A*x - b), both computed from the explicit
    %   matrices on a small mesh. The preconditioned inverse iteration pinvit with a
    %   multigrid preconditioner finds the smallest eigenvalue of the
    %   discrete Laplacian with a B-normalized eigenvector.
    %
    %   See also FE, FEVector, Multigrid, pinvit, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degree = struct('P1',1,'P2',2,'P3',3)
        degree3d = struct('P1',1,'P2',2)
        refinement = struct('uniform','uniform','RGB','RGB','NVB','NVB')
        element2d = struct( ...
            'CR1',struct('space',@(m)CRk(1,m),'field','scalar','fieldDegree',1), ...
            'CR3',struct('space',@(m)CRk(3,m),'field','scalar','fieldDegree',3), ...
            'DGP2',struct('space',@(m)DGPk(2,m),'field','scalar','fieldDegree',2), ...
            'DGD2',struct('space',@(m)DGDk(2,m),'field','scalar','fieldDegree',2), ...
            'HP3',struct('space',@(m)HPk(3,m),'field','scalar','fieldDegree',3), ...
            'RT1',struct('space',@(m)RTk(1,m),'field','vector','fieldDegree',0), ...
            'RT3',struct('space',@(m)RTk(3,m),'field','vector','fieldDegree',2), ...
            'BDM2',struct('space',@(m)BDMk(2,m),'field','vector','fieldDegree',2), ...
            'Nd1_2',struct('space',@(m)Nd1k(2,m),'field','vector','fieldDegree',1), ...
            'Nd2_2',struct('space',@(m)Nd2k(2,m),'field','vector','fieldDegree',2), ...
            'HHJ1',struct('space',@(m)HHJk(1,m),'field','symmatrix','fieldDegree',1), ...
            'TUBA5',struct('space',@(m)TUBAk(5,m),'field','scalar','fieldDegree',5), ...
            'Argyris',struct('space',@(m)Argyris(5,m),'field','scalar','fieldDegree',5), ...
            'AW1',struct('space',@(m)AW1(3,m),'field','symmatrix','fieldDegree',2), ...
            'HZ1',struct('space',@(m)HZ1(3,m),'field','symmatrix','fieldDegree',2), ...
            'AWNC1',struct('space',@(m)AWNC1(2,m),'field','symmatrix','fieldDegree',1), ...
            'RBM',struct('space',@(m)RBM(1,m),'field','rigid','fieldDegree',0))
        refinement2d = struct( ...
            'RGB',@(m)m.refine(MultilevelTest.everyThird(m),'RGB'), ...
            'NVB',@(m)m.refine(MultilevelTest.everyThird(m),'NVB'), ...
            'uniform2',@(m)m.refine('uniform',2), ...
            'graded',@(m)m.refine('graded',0.5,0.3,[0 0]), ...
            'gradedRed',@(m)m.refine('gradedred',0.3,[0 0]))
        orientedElement = struct( ...
            'RT1',struct('name','RTk3d','degree',1,'field','vector','fieldDegree',0), ...
            'RT2',struct('name','RTk3d','degree',2,'field','vector','fieldDegree',1), ...
            'BDM1',struct('name','BDMk3d','degree',1,'field','vector','fieldDegree',1), ...
            'BDM2',struct('name','BDMk3d','degree',2,'field','vector','fieldDegree',2), ...
            'Nd2_1',struct('name','Nd2k3d','degree',1,'field','vector','fieldDegree',1), ...
            'Nd2_2',struct('name','Nd2k3d','degree',2,'field','vector','fieldDegree',2), ...
            'Nd1_1',struct('name','Nd1k3d','degree',1,'field','vector','fieldDegree',0), ...
            'Nd1_2',struct('name','Nd1k3d','degree',2,'field','vector','fieldDegree',1), ...
            'HP2',struct('name','HPk3d','degree',2,'field','scalar','fieldDegree',2))
        marking = struct('partial','partial','graded','graded','uniform','uniform')
        sumSpace = struct( ...
            'BDM1CurlBRot2',struct('space',@(m)FEMatrix(BDMk(1,m),BDMk(1,m)) + CurlBRotk(2,m), ...
            'first',@(m)FEMatrix(BDMk(1,m),BDMk(1,m)),'nested',true), ...
            'RT2CurlBRot2',struct('space',@(m)FEMatrix(RTk(2,m),RTk(2,m)) + CurlBRotk(2,m), ...
            'first',@(m)FEMatrix(RTk(2,m),RTk(2,m)),'nested',true), ...
            'BDM2CurlBRot3',struct('space',@(m)FEMatrix(BDMk(2,m),BDMk(2,m)) + CurlBRotk(3,m), ...
            'first',@(m)FEMatrix(BDMk(2,m),BDMk(2,m)),'nested',true), ...
            'Stenberg1',struct('space',@(m)FEMatrix(BDMk(1,m),BDMk(1,m)) + FEMatrix(CurlBk(2,m),CurlBk(2,m)), ...
            'first',@(m)FEMatrix(BDMk(1,m),BDMk(1,m)),'nested',true), ...
            'Mini',struct('space',@(m)PkVector(1,m) + FEVector(Bk(3,m),Bk(3,m)), ...
            'first',@(m)PkVector(1,m),'nested',false), ...
            'Mini3to4',struct('space',@(m)PkVector(1,m) + FEVector(Bk(3+mod((1:m.nrElems)',2),m), ...
            Bk(3+mod((1:m.nrElems)',2),m)),'first',@(m)PkVector(1,m),'nested',false))
        sameSpace = struct( ...
            'P3',struct('space',@(m)Pk(3,m),'dim',2,'affine',true), ...
            'CR3',struct('space',@(m)CRk(3,m),'dim',2,'affine',true), ...
            'DGP2',struct('space',@(m)DGPk(2,m),'dim',2,'affine',true), ...
            'DGD2',struct('space',@(m)DGDk(2,m),'dim',2,'affine',true), ...
            'HP3',struct('space',@(m)HPk(3,m),'dim',2,'affine',true), ...
            'DGHP3',struct('space',@(m)DGHPk(3,m),'dim',2,'affine',true), ...
            'RT2',struct('space',@(m)RTk(2,m),'dim',2,'affine',true), ...
            'HRT2',struct('space',@(m)HRTk(2,m),'dim',2,'affine',true), ...
            'BDM2',struct('space',@(m)BDMk(2,m),'dim',2,'affine',true), ...
            'HBDM1',struct('space',@(m)HBDMk(1,m),'dim',2,'affine',true), ...
            'Nd1_2',struct('space',@(m)Nd1k(2,m),'dim',2,'affine',true), ...
            'Nd2_2',struct('space',@(m)Nd2k(2,m),'dim',2,'affine',true), ...
            'HHJ1',struct('space',@(m)HHJk(1,m),'dim',2,'affine',true), ...
            'HHJB1',struct('space',@(m)HHJB(1,m),'dim',2,'affine',true), ...
            'B4',struct('space',@(m)Bk(4,m),'dim',2,'affine',true), ...
            'B2_6',struct('space',@(m)B2k(6,m),'dim',2,'affine',true), ...
            'CurlB3',struct('space',@(m)CurlBk(3,m),'dim',2,'affine',true), ...
            'CurlBRot2',struct('space',@(m)CurlBRotk(2,m),'dim',2,'affine',true), ...
            'TUBA5',struct('space',@(m)TUBAk(5,m),'dim',2,'affine',false), ...
            'Argyris',struct('space',@(m)Argyris(5,m),'dim',2,'affine',false), ...
            'AW1',struct('space',@(m)AW1(3,m),'dim',2,'affine',false), ...
            'AW0',struct('space',@(m)AW0(3,m),'dim',2,'affine',false), ...
            'HZ1',struct('space',@(m)HZ1(3,m),'dim',2,'affine',false), ...
            'AW4',struct('space',@(m)AWk(4,m),'dim',2,'affine',false), ...
            'AWH3',struct('space',@(m)AWH(3,m),'dim',2,'affine',false), ...
            'AWNC1',struct('space',@(m)AWNC1(2,m),'dim',2,'affine',false), ...
            'AWNC0',struct('space',@(m)AWNC0(2,m),'dim',2,'affine',false), ...
            'GGNC1',struct('space',@(m)GGNC1(2,m),'dim',2,'affine',false), ...
            'RBM',struct('space',@(m)RBM(1,m),'dim',2,'affine',false), ...
            'PVector2',struct('space',@(m)PkVector(2,m),'dim',2,'affine',true), ...
            'DGPVector1',struct('space',@(m)DGPkVector(1,m),'dim',2,'affine',true), ...
            'DGDMatrix1',struct('space',@(m)DGDkMatrix(1,m),'dim',2,'affine',true), ...
            'DGPSymMatrix1',struct('space',@(m)DGPkSymMatrix(1,m),'dim',2,'affine',true), ...
            'P2_3d',struct('space',@(m)Pk3d(2,m),'dim',3,'affine',true), ...
            'DGP2_3d',struct('space',@(m)DGPk3d(2,m),'dim',3,'affine',true), ...
            'DGD2_3d',struct('space',@(m)DGDk3d(2,m),'dim',3,'affine',true), ...
            'HP2_3d',struct('space',@(m)HPk3d(2,m),'dim',3,'affine',true), ...
            'RT2_3d',struct('space',@(m)RTk3d(2,m),'dim',3,'affine',true), ...
            'BDM1_3d',struct('space',@(m)BDMk3d(1,m),'dim',3,'affine',true), ...
            'Nd1_2_3d',struct('space',@(m)Nd1k3d(2,m),'dim',3,'affine',true), ...
            'Nd2_1_3d',struct('space',@(m)Nd2k3d(1,m),'dim',3,'affine',true), ...
            'Nd1DG2_3d',struct('space',@(m)Nd1DGk3d(2,m),'dim',3,'affine',true), ...
            'PVector1_3d',struct('space',@(m)PkVector(1,m),'dim',3,'affine',true))
        affineElement = struct( ...
            'P3',@(m)Pk(3,m),'CR3',@(m)CRk(3,m),'DGD2',@(m)DGDk(2,m),'HP3',@(m)HPk(3,m), ...
            'RT2',@(m)RTk(2,m),'BDM2',@(m)BDMk(2,m),'Nd1_2',@(m)Nd1k(2,m),'Nd2_2',@(m)Nd2k(2,m), ...
            'HHJ1',@(m)HHJk(1,m),'B4',@(m)Bk(4,m))
    end

    methods (Test)

        function testProlongation(testCase,degree,refinement)
            coarse = Mesh(struct('name','Square')).refine('uniform');
            if strcmp(refinement,'uniform')
                fine = coarse.refine('uniform');
            else
                fine = coarse.refine(MultilevelTest.everyThird(coarse),refinement);
            end

            MultilevelTest.verifyProlongation(testCase,coarse,fine,@(m)Pk(degree,m), ...
                @(fe,feOld)fe.prolongation(feOld));
        end

        function testGradedRedProlongation(testCase,degree)
            coarse = Mesh(struct('name','Square')).refine('uniform');

            fine = coarse.refine('gradedred',0.3,[0 0]);

            MultilevelTest.verifyProlongation(testCase,coarse,fine,@(m)Pk(degree,m), ...
                @(fe,feOld)fe.prolongation(feOld));
        end

        function testElementProlongation2d(testCase,element2d,refinement2d)
            coarse = Mesh(struct('name','Square')).refine('uniform');
            fine = refinement2d(coarse);
            feOld = element2d.space(coarse);
            fe = element2d.space(fine);
            f = ProjectionTest.field(element2d.field,element2d.fieldDegree,2);
            x = fe.project(f);

            P = fe.prolongation(feOld);

            testCase.verifySize(P,[fe.nrDoF,feOld.nrDoF]);
            testCase.verifyEqual(P*feOld.project(f),x,'AbsTol',1e-8*max(1,norm(x,inf)));
        end

        function testEnrichedProlongation(testCase)
            % HHJB, HHJk with the TDNNS bubbles, contains HHJk of the same
            % degree
            coarse = Mesh(struct('name','Square')).refine('uniform');
            fine = coarse.refine(MultilevelTest.everyThird(coarse),'RGB');
            feOld = HHJk(1,coarse);
            fe = HHJB(1,fine);
            f = ProjectionTest.field('symmatrix',1,2);

            P = fe.prolongation(feOld);

            testCase.verifySize(P,[fe.nrDoF,feOld.nrDoF]);
            testCase.verifyEqual(P*feOld.project(f),fe.project(f),'AbsTol',1e-8);
        end

        function testBubbleProlongation(testCase)
            % bubble spaces, alone and in the Mini element, and HHJB, which
            % the reference prolongation rebuilds by (degree,mesh)
            coarse = Mesh(struct('name','Square')).refine('uniform');
            fine = coarse.refine(MultilevelTest.everyThird(coarse),'RGB');
            spaces = {@(m)Bk(3,m),@(m)Bk(4,m),@(m)B2k(6,m),@(m)CurlBk(3,m),@(m)CurlBRotk(2,m), ...
                @(m)PkVector(1,m)+FEVector(Bk(3,m),Bk(3,m)),@(m)HHJB(2,m)};
            for k = 1 : numel(spaces)
                fe = spaces{k}(fine);
                feOld = spaces{k}(coarse);

                P = fe.prolongation(feOld);

                testCase.verifySize(P,[fe.nrDoF,feOld.nrDoF],func2str(spaces{k}));
                testCase.verifyTrue(all(isfinite(nonzeros(P))),func2str(spaces{k}));
            end
        end

        function testSumProlongation(testCase,sumSpace,refinement)
            % local L2 projection onto the whole fine sum space
            coarse = Mesh(struct('name','Square')).refine('uniform');
            if strcmp(refinement,'uniform')
                fine = coarse.refine('uniform');
            else
                fine = coarse.refine(MultilevelTest.everyThird(coarse),refinement);
            end
            feOld = sumSpace.space(coarse);
            fe = sumSpace.space(fine);
            rng(3);
            xOld = randn(feOld.allnrDoF,1);
            xOld(feOld.nonDoFNodes) = 0;
            if ~sumSpace.nested || ~strcmp(refinement,'uniform')
                % only the first summand is nested
                xOld(sumSpace.first(coarse).allnrDoF+1:end) = 0;
            end

            P = fe.prolongation(feOld);

            testCase.verifySize(P,[fe.allnrDoF,feOld.allnrDoF]);
            unused = setdiff(1:fe.allnrDoF,[fe.freeNodes(:);fe.fixedNodes(:)]);
            testCase.verifyEqual(nnz(P(unused,:)),0,'unused and linearly dependent dofs');
            % the fine function equals the coarse one
            ProjectionTest.verifyValues(testCase,fe,P*xOld, ...
                @(z)ProjectionTest.valueInParent(feOld,xOld,fine.point4e(z),fine.parent4e));
        end

        function testCurlBRotkMismatch(testCase)
            % FEMatrix + CurlBRotk of the weakly symmetric elasticity needs
            % a coarse space with a CurlBRotk element, and vice versa
            coarse = Mesh(struct('name','Square')).refine('uniform');
            fine = coarse.refine('uniform');
            space = @(m)FEMatrix(BDMk(1,m),BDMk(1,m));

            testCase.verifyError(@()prolongation(space(fine) + CurlBRotk(2,fine),space(coarse)), ...
                'fem4m:FEMatrix:curlBRotkMismatch');
            testCase.verifyError(@()space(fine).prolongation(space(coarse) + CurlBRotk(2,coarse)), ...
                'fem4m:FEMatrix:curlBRotkMismatch');
        end

        function testElementwiseVectorProlongation(testCase,refinement)
            % components with elementwise degrees 0 to 2, inherited from
            % the parents
            coarse = Mesh(struct('name','Square')).refine('uniform');
            if strcmp(refinement,'uniform')
                fine = coarse.refine('uniform');
            else
                fine = coarse.refine(MultilevelTest.everyThird(coarse),refinement);
            end
            degree4e = mod((1:coarse.nrElems)',3);
            feOld = DGDkVector(degree4e,coarse);
            fe = DGDkVector(degree4e(fine.parent4e),fine);
            rng(5);
            xOld = randn(feOld.allnrDoF,1);
            xOld(feOld.nonDoFNodes) = 0;

            P = fe.prolongation(feOld);

            testCase.verifySize(P,[fe.allnrDoF,feOld.allnrDoF]);
            testCase.verifyEqual(nnz(P(fe.nonDoFNodes,:)),0,'unused dofs');
            ProjectionTest.verifyValues(testCase,fe,P*xOld, ...
                @(z)ProjectionTest.valueInParent(feOld,xOld,fine.point4e(z),fine.parent4e));
        end

        function testRebuildFromDegree(testCase,sameSpace)
            % the reference prolongation builds an element on other
            % meshes by feval(class(fe),fe.degree,mesh)
            if sameSpace.dim == 2
                mesh = Mesh(struct('name','Square')).refine('uniform');
                ref = [0.2 0.3];
            else
                mesh = UnitTestCase.cubeMesh();
                ref = [0.2 0.3 0.1];
            end
            fe = sameSpace.space(mesh);
            parts = 1:mesh.nrElems;

            other = feval(class(fe),fe.degree,mesh);

            testCase.verifyClass(other,class(fe));
            testCase.verifyEqual([other.nrBasisFunc,other.nrDoF],[fe.nrBasisFunc,fe.nrDoF]);
            testCase.verifyEqual(other.dof4e,fe.dof4e);
            testCase.verifyEqual(other.getBasis(mesh.point4e(ref),ref,parts), ...
                fe.getBasis(mesh.point4e(ref),ref,parts));
            testCase.verifyEqual(fe.isAffineEquivalent,sameSpace.affine,'isAffineEquivalent');
        end

        function testFactoredEqualsPhysical(testCase,affineElement)
            % for affine equivalent elements the reference matrices give
            % the local L2 projections on the physical elements, except
            % for Piola-mapped spaces that are not nested (CurlBk,
            % CurlBRotk, HHJB), where they project in reference
            % coordinates
            coarse = Mesh(struct('name','Square')).refine('uniform');
            fine = coarse.refine(MultilevelTest.everyThird(coarse),'RGB');
            feOld = affineElement(coarse);
            fe = affineElement(fine);
            P = fe.prolongation(feOld);
            fe.isAffineEquivalent = false;

            Pphysical = fe.prolongation(feOld);

            testCase.verifyEqual(full(Pphysical),full(P),'AbsTol',1e-10*max(1,max(abs(P(:)))));
        end

        function testSmallScaleProlongation(testCase)
            % on elements of size 1e-6 the local matrices of Argyris
            % couple derivative dofs by entries of size h^2
            s = 1e-6;
            unit = Mesh(struct('name','Square')).refine('uniform');
            coarse = Mesh(s*unit.c4n,unit.n4e,unit.n4Db,unit.n4Nb);
            fine = coarse.refine(MultilevelTest.everyThird(coarse),'RGB');
            p = UnitTestCase.polynomial(5);
            feOld = Argyris(5,coarse);
            fe = Argyris(5,fine);
            % no warnings, although the unscaled local mass matrices of
            % project have rcond 1e-33
            xOld = testCase.verifyWarningFree(@()feOld.project(@(x)p(x/s)));
            P = testCase.verifyWarningFree(@()fe.prolongation(feOld));

            x = P*xOld;

            % the fine function equals the coarse one at the centroids
            centroid = fine.point4e([1/3,1/3]);
            value = fe.getValue([],[1/3,1/3],[],x);
            coarseValue = zeros(fine.nrElems,1);
            for e = 1 : fine.nrElems
                parent = fine.parent4e(e);
                C = coarse.c4n(coarse.n4e(parent,:),:);
                ref = ([C(2,:)-C(1,:);C(3,:)-C(1,:)]'\(centroid(e,:)-C(1,:))')';
                coarseValue(e) = feOld.getValue([],ref,parent,xOld);
            end
            testCase.verifyEqual(value,coarseValue,'AbsTol',1e-8*max(1,norm(coarseValue,inf)));
        end

        function testHighDegreeProlongation(testCase)
            % element-wise degrees up to 12 (local L2 projections with the
            % coarse basis fitted in the Dubiner basis)
            coarse = Mesh(struct('name','Square')).refine('uniform');
            fine = coarse.refine(MultilevelTest.everyThird(coarse),'RGB');
            degree4e = 12 - mod((1:coarse.nrElems)',3);
            feOld = DGDk(degree4e,coarse);
            fe = DGDk(degree4e(fine.parent4e),fine);
            f = UnitTestCase.polynomial(10);
            x = fe.project(f);

            P = fe.prolongation(feOld);

            testCase.verifyEqual(P*feOld.project(f),x,'AbsTol',1e-8*max(1,norm(x,inf)));
        end

        function testStronglyGradedProlongation(testCase)
            % fine elements far smaller than their parents in one refine
            % call: the local L2 projections keep their small entries, and
            % the fine function equals the coarse one
            coarse = Mesh(struct('name','Square')).refine('uniform');
            fine = coarse.refine('graded',0.2,0.3,[0 0]);
            degree4e = 8 - mod((1:coarse.nrElems)',3);
            feOld = DGDk(degree4e,coarse);
            fe = DGDk(degree4e(fine.parent4e),fine);
            rng(7);
            xOld = randn(feOld.allnrDoF,1);
            xOld(feOld.nonDoFNodes) = 0;

            x = fe.prolongation(feOld)*xOld;

            centroid = fine.point4e([1/3,1/3]);
            value = fe.getValue([],[1/3,1/3],[],x);
            coarseValue = zeros(fine.nrElems,1);
            for e = 1 : fine.nrElems
                parent = fine.parent4e(e);
                C = coarse.c4n(coarse.n4e(parent,:),:);
                ref = ([C(2,:)-C(1,:);C(3,:)-C(1,:)]'\(centroid(e,:)-C(1,:))')';
                coarseValue(e) = feOld.getValue([],ref,parent,xOld);
            end
            testCase.verifyEqual(value,coarseValue,'AbsTol',1e-12*norm(coarseValue,inf));
        end

        function testIdentityProlongation(testCase)
            % no refinement round: the local L2 projections give the
            % identity without round-off entries
            mesh = Mesh(struct('name','Square')).refine('uniform');
            same = mesh.refine('uniform',0);
            spaces = {@(m)Argyris(5,m),@(m)TUBAk(8,m),@(m)AW1(3,m)};
            for k = 1 : numel(spaces)
                fe = spaces{k}(same);

                P = fe.prolongation(spaces{k}(mesh));

                testCase.verifyEqual(nnz(P),fe.nrDoF,func2str(spaces{k}));
                testCase.verifyEqual(full(diag(P)),ones(fe.nrDoF,1),'AbsTol',1e-10,func2str(spaces{k}));
            end
        end

        function testSkewSymProlongation(testCase,refinement)
            % the basis of FESkewSymMatrix is scaled by the patch areas
            coarse = Mesh(struct('name','Square')).refine('uniform');
            if strcmp(refinement,'uniform')
                fine = coarse.refine('uniform');
            else
                fine = coarse.refine(MultilevelTest.everyThird(coarse),refinement);
            end
            p = UnitTestCase.polynomial(1);
            f = @(z)[0*z(:,1),p(z),-p(z),0*z(:,1)];
            spaces = {@(m)DGDk(1,m),@(m)Pk(1,m)};
            for k = 1 : numel(spaces)
                feOld = FESkewSymMatrix(spaces{k}(coarse));
                fe = FESkewSymMatrix(spaces{k}(fine));
                x = fe.project(f);

                P = fe.prolongation(feOld);

                testCase.verifyEqual(P*feOld.project(f),x,'AbsTol',1e-8*max(1,norm(x,inf)),func2str(spaces{k}));
                testCase.verifyEqual(fe.prolongation(spaces{k}(coarse)),P,'scalar element');
            end
        end

        function testRefineParents2d(testCase)
            coarse = Mesh(struct('name','Square')).refine('uniform');

            fine = coarse.refine('uniform',2);

            % parent4e and the P1 prolongation refer to the mesh before
            % the call: the fine nodes are the P1 interpolants of the
            % coordinates, and every element lies in its ancestor
            testCase.verifySize(fine.prolongation,[fine.nrNodes,coarse.nrNodes]);
            testCase.verifyEqual(fine.prolongation*coarse.c4n,fine.c4n,'AbsTol',testCase.tol);
            testCase.verifyEqual(numel(fine.parent4e),fine.nrElems);
            testCase.verifyEqual(accumarray(fine.parent4e(:),fine.area4e), ...
                coarse.area4e,'AbsTol',testCase.tol);
        end

        function testDegreeIncrease(testCase)
            coarse = Mesh(struct('name','Square')).refine('uniform');
            fine = coarse.refine(MultilevelTest.everyThird(coarse),'NVB');
            feOld = Pk(1,coarse);
            fe = Pk(2,fine);
            f = MultilevelTest.polynomial(1,2);

            P = fe.prolongation(feOld);

            testCase.verifyEqual(P*feOld.project(f),fe.project(f),'AbsTol',testCase.tol);
        end

        function testTwoRefineCalls(testCase)
            mesh0 = Mesh(struct('name','Lshape'));
            mesh1 = mesh0.refine(MultilevelTest.everyThird(mesh0),'NVB');
            mesh2 = mesh1.refine(MultilevelTest.everyThird(mesh1),'RGB');
            fe0 = RTk(2,mesh0);
            fe1 = RTk(2,mesh1);
            fe2 = RTk(2,mesh2);
            f = ProjectionTest.field('vector',1,2);

            P = fe2.prolongation(fe1)*fe1.prolongation(fe0);

            testCase.verifyEqual(P*fe0.project(f),fe2.project(f),'AbsTol',testCase.tol);
            testCase.verifyError(@()fe2.prolongation(fe0),'fem4m:FE:notParentMesh');
        end

        function testVectorProlongation2d(testCase)
            coarse = Mesh(struct('name','Square')).refine('uniform');
            fine = coarse.refine(MultilevelTest.everyThird(coarse),'RGB');
            feOld = PkVector(2,coarse);
            fe = PkVector(2,fine);
            f = ProjectionTest.field('vector',2,2);
            xOld = feOld.project(f);

            P = fe.prolongation(feOld);

            testCase.verifySize(P,[fe.nrDoF,feOld.nrDoF]);
            testCase.verifyEqual(P*xOld,fe.project(f),'AbsTol',testCase.tol);
            testCase.verifyEqual(fe.prolongate(xOld,Pk(2,coarse)),P*xOld, ...
                'AbsTol',testCase.tol,'one scalar coarse element for all components');
        end

        function testProlongation3d(testCase,degree3d,marking)
            [coarse,fine] = MultilevelTest.meshes3d(marking);

            MultilevelTest.verifyProlongation(testCase,coarse,fine,@(m)Pk3d(degree3d,m), ...
                @(fe,feOld)fe.prolongation(feOld));
        end

        function testOrientedProlongation3d(testCase,orientedElement,marking)
            space = @(m)feval(orientedElement.name,orientedElement.degree,m);
            [coarse,fine] = MultilevelTest.meshes3d(marking);
            feOld = space(coarse);
            fe = space(fine);
            f = ProjectionTest.field(orientedElement.field,orientedElement.fieldDegree,3);

            P = fe.prolongation(feOld);

            testCase.verifySize(P,[fe.nrDoF,feOld.nrDoF]);
            testCase.verifyEqual(P*feOld.project(f),fe.project(f),'AbsTol',testCase.tol);
        end

        function testContractionNumber(testCase)
            [mg,A] = MultilevelTest.poissonMultigrid(3);
            n = size(A,1);

            rho = mg.contractionNumber;

            % explicit error propagation of one V-cycle from zero
            BA = zeros(n);
            for j = 1 : n
                BA(:,j) = mg.solve(A(:,j),0,1,[],0);
            end
            expected = max(abs(eig(eye(n) - BA)));
            testCase.verifyGreaterThan(expected,0);
            testCase.verifyLessThan(expected,1,'convergent multigrid');
            testCase.verifyEqual(rho,expected,'RelTol',1e-2);
        end

        function testDampingFactor(testCase)
            % with the damping 1/16, the spectra of T*A and I - T*A
            % differ
            [mg,A] = MultilevelTest.poissonMultigrid(3,1/16);
            n = size(A,1);
            damp = mg.damp{end};

            val = mg.dampingFactor;

            % Richardson smoother x - damp*(A*x - b), i.e. T = damp*I
            expected = max(abs(eig(damp.*full(A))));
            testCase.assertGreaterThan(abs(expected - max(abs(eig(eye(n) - damp.*full(A))))),0.1);
            testCase.verifyEqual(val,expected,'RelTol',1e-2);
        end

        function testPinvit(testCase)
            [mg,A,B] = MultilevelTest.poissonMultigrid(4);
            x0 = ones(size(A,1),1);

            [x,lambda,iter] = pinvit(x0,A,B,1e-10,mg,2,200);

            expected = eigs(A,B,1,'smallestabs');
            testCase.verifyLessThan(iter,200,'converged');
            testCase.verifyEqual(lambda,expected,'RelTol',1e-8,'smallest eigenvalue');
            testCase.verifyEqual(x'*B*x,1,'AbsTol',testCase.tol,'B-normalized');
            testCase.verifyLessThan(norm(A*x - lambda*B*x)/norm(B*x),1e-4,'eigenvector');
        end

    end

    methods (Static)

        function verifyProlongation(testCase,coarse,fine,space,prolongation)
            %verifyProlongation - Prolongation of a polynomial of degree k
            feOld = space(coarse);
            fe = space(fine);
            f = MultilevelTest.polynomial(fe.degree,size(coarse.c4n,2));
            xOld = feOld.project(f);
            x = fe.project(f);

            P = prolongation(fe,feOld);

            testCase.verifySize(P,[fe.nrDoF,feOld.nrDoF]);
            testCase.verifyEqual(P*xOld,x,'AbsTol',testCase.tol);
        end

        function [coarse,fine] = meshes3d(marking)
            %meshes3d - Coarse and refined cube meshes for 3d prolongations
            %   'partial' bisects one element of a mesh with elements of
            %   two generations, so that 12 of its 14 elements stay
            %   unrefined; 'graded' starts from that mesh and bisects all
            %   elements, some of them twice in one call; 'uniform'
            %   refines the cube uniformly.
            switch marking
                case 'partial'
                    coarse = UnitTestCase.cubeMesh();
                    fine = coarse.refine([true;false(coarse.nrElems-1,1)]);
                case 'graded'
                    coarse = UnitTestCase.cubeMesh();
                    fine = coarse.refine(true(coarse.nrElems,1));
                otherwise
                    coarse = Mesh3d(Geometry(struct('name','Cube')));
                    fine = coarse.refine('uniform');
            end
        end

        function marked4e = everyThird(mesh)
            %everyThird - Mark every third element
            marked4e = zeros(mesh.nrElems,1);
            marked4e(1:3:end) = 1;
        end

        function f = polynomial(k,dim)
            %polynomial - Polynomial of total degree k in 2d or 3d
            if dim == 2
                f = UnitTestCase.polynomial(k);
            else
                f = @(z)z(:,1).^k - 2*z(:,1).*z(:,2).^(k-1) + z(:,3).^k + z(:,2) + 1;
            end
        end

        function [mg,A,B] = poissonMultigrid(nrLevels,damp)
            %poissonMultigrid - P1 multigrid for the Laplacian on the square
            %   A and B are the stiffness and mass matrices of the free
            %   nodes on the finest level. DAMP is the damping of the
            %   Richardson smoother, 1/8 by default.
            if nargin < 2
                damp = [];
            end
            mesh = Mesh(struct('name','Square'));
            mg = Multigrid('Richardson',damp);
            pk = [];
            for k = 1 : nrLevels
                mesh = mesh.refine('uniform');
                pkOld = pk;
                pk = Pk(1,mesh);
                A = BilinearForm('grad(u) * grad(v) dx').assemble(pk);
                A = A(pk.freeNodes,pk.freeNodes);
                if k == 1
                    mg = mg.addLevel(A);
                else
                    P = pk.prolongation(pkOld);
                    mg = mg.addLevel(A,P(pk.freeNodes,pkOld.freeNodes));
                end
            end
            B = BilinearForm('u * v dx').assemble(pk);
            B = B(pk.freeNodes,pk.freeNodes);
        end

    end

end
classdef BubbleElementTest < UnitTestCase
    %BubbleElementTest - Unit tests of the curls of bubble functions
    %   The basis of CurlBk(k) is the curl (d_y phi, -d_x phi) of the
    %   basis phi of Bk(k+1), the cubic bubble b_T times the polynomials
    %   of degree k-2, with the gradients of the curls given by the
    %   Hessians of phi and vanishing divergence. The basis of
    %   CurlBRotk(k) holds the matrices whose rows are the curls
    %   (d_y psi, -d_x psi) of psi = b_T d_x q and psi = b_T d_y q for
    %   the non-constant Dubiner basis functions q of degree k-1 >= 1,
    %   scaled by the square root of the element area, so that the row
    %   divergence vanishes. HHJB(k) adds 2k divdiv-free bubbles to the
    %   basis of HHJk(k). The basis functions of a finite element are
    %   linearly independent.
    %
    %   See also CurlBk, CurlBRotk, Bk, DGDk, HHJB, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degree = struct('k2',2,'k3',3,'k4',4)
        rotDegree = struct('k2',2,'k3',3,'k4',4)
        hhjDegree = struct('k1',1,'k2',2)
    end

    methods (Test)

        function testCurlBubble(testCase,degree)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = CurlBk(degree,mesh);
            bubble = Bk(degree+1,mesh);
            ref = [0.2 0.3];
            pts = mesh.point4e(ref);
            parts = 1:mesh.nrElems;

            val = BubbleElementTest.rows(fe,'getBasis',pts,ref);
            grad = BubbleElementTest.rows(fe,'getGradBasis',pts,ref);
            div = BubbleElementTest.rows(fe,'getDivBasis',pts,ref);

            g = BubbleElementTest.rows(bubble,'getGradBasis',pts,ref);
            H = BubbleElementTest.rows(bubble,'getD2Basis',pts,ref);  % [xx yy xy]
            testCase.verifyEqual(fe.nrBasisFunc,bubble.nrBasisFunc);
            testCase.verifyEqual(fe.degree,degree);
            testCase.verifyEqual(val,cat(3,g(:,:,2),-g(:,:,1)),'AbsTol',testCase.tol,'curl');
            % J(i,j) = d_j u_i column by column
            testCase.verifyEqual(grad,cat(3,H(:,:,3),-H(:,:,1),H(:,:,2),-H(:,:,3)), ...
                'AbsTol',testCase.tol,'gradient of the curl');
            testCase.verifyEqual(grad(:,:,1) + grad(:,:,4),zeros(numel(parts),fe.nrBasisFunc), ...
                'AbsTol',testCase.tol,'divergence free');
            testCase.verifyEqual(div,zeros(numel(parts),fe.nrBasisFunc));
        end

        function testCurlBubbleRotation(testCase,rotDegree)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = CurlBRotk(rotDegree,mesh);
            dg = DGDk(rotDegree-1,mesh);
            cubic = Bk(3,mesh);
            ref = [0.2 0.3];
            pts = mesh.point4e(ref);

            val = BubbleElementTest.rows(fe,'getBasis',pts,ref);

            % b_T = cubic/cubic(1/3,1/3) is the only basis function of Bk(3)
            scale = BubbleElementTest.rows(cubic,'getBasis',mesh.point4e([1 1]/3),[1 1]/3)*27;
            b = BubbleElementTest.rows(cubic,'getBasis',pts,ref)./scale;
            gb = BubbleElementTest.rows(cubic,'getGradBasis',pts,ref)./scale;
            gq = BubbleElementTest.rows(dg,'getGradBasis',pts,ref);
            Hq = BubbleElementTest.rows(dg,'getD2Basis',pts,ref);  % [xx yy xy]
            % curls of b d_x q (first row) and b d_y q (second row)
            M11 = gb(:,:,2).*gq(:,:,1) + b.*Hq(:,:,3);
            M12 = -gb(:,:,1).*gq(:,:,1) - b.*Hq(:,:,1);
            M21 = gb(:,:,2).*gq(:,:,2) + b.*Hq(:,:,2);
            M22 = -gb(:,:,1).*gq(:,:,2) - b.*Hq(:,:,3);
            % without the vanishing function of the constant polynomial
            expected = sqrt(mesh.area4e).*cat(3,M11,M21,M12,M22);
            expected = expected(:,2:end,:);
            testCase.verifyEqual(fe.nrBasisFunc,dg.nrBasisFunc - 1);
            testCase.verifyEqual(val,expected,'AbsTol',testCase.tol);
            testCase.verifyEqual(BubbleElementTest.rows(fe,'getDivBasis',pts,ref), ...
                zeros(mesh.nrElems,fe.nrBasisFunc,2));
        end

        function testDivDivFreeBubbles(testCase,hhjDegree)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = HHJB(hhjDegree,mesh);
            hhj = HHJk(hhjDegree,mesh);
            ref = [0.2 0.3];
            pts = mesh.point4e(ref);

            basis = BubbleElementTest.rows(fe,'getBasis',pts,ref);
            divdiv = BubbleElementTest.rows(fe,'getDivDivBasis',pts,ref);

            testCase.verifyTrue(fe.bubbles);
            testCase.verifyEqual(fe.nrBasisFunc,hhj.nrBasisFunc + 2*hhjDegree);
            testCase.verifyEqual(basis(:,1:hhj.nrBasisFunc,:), ...
                BubbleElementTest.rows(hhj,'getBasis',pts,ref),'AbsTol',testCase.tol, ...
                'HHJk basis first');
            testCase.verifyEqual(divdiv(:,hhj.nrBasisFunc+1:end), ...
                zeros(mesh.nrElems,2*hhjDegree),'AbsTol',testCase.tol,'divdiv free');
            testCase.verifyEqual(BubbleElementTest.rank(fe),fe.nrBasisFunc);
        end

        function testLinearIndependence(testCase,degree)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));

            fe = CurlBk(degree,mesh);

            testCase.verifyEqual(BubbleElementTest.rank(fe),fe.nrBasisFunc);
        end

        function testLinearIndependenceRotation(testCase,rotDegree)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));

            fe = CurlBRotk(rotDegree,mesh);

            testCase.verifyEqual(BubbleElementTest.rank(fe),fe.nrBasisFunc);
        end

    end

    methods (Static)

        function r = rank(fe)
            %rank - Rank of the basis of the first element at many points
            [x,y] = meshgrid(linspace(0.05,0.9,7));
            inside = x + y < 0.95;
            refPoints = [x(inside),y(inside)];
            V = [];
            for k = 1 : size(refPoints,1)
                ref = refPoints(k,:);
                basis = BubbleElementTest.rows(fe,'getBasis',fe.mesh.point4e(ref),ref);
                V = [V;reshape(basis(1,:,:),fe.nrBasisFunc,[])']; %#ok<AGROW>
            end
            r = rank(V);
        end

        function val = rows(fe,method,pts,ref)
            %rows - Basis METHOD of all elements at one reference point
            parts = 1:fe.mesh.nrElems;
            val = FormBasis.elementRows(fe.(method)(pts,ref,parts),parts,fe);
        end

    end

end
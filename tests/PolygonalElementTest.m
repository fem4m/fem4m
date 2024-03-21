classdef PolygonalElementTest < UnitTestCase
    %PolygonalElementTest - Unit tests of the scaled monomial spaces
    %   PDGk and its vector, matrix and symmetric matrix versions are
    %   spanned by the monomials scaled by the centroid and diameter of the
    %   polygons of MeshPolyT, MeshPolyQ and MeshPolyH. The L2 projection
    %   reproduces the polynomial fields they contain and their gradients,
    %   also with a different scaling on every element, and the row
    %   divergence, row-wise rot and trace of PDGkSymMatrix are those of the
    %   symmetric matrix field. The second, third and fourth derivative
    %   bases of PDGk give those of a sextic polynomial, and the Hessian,
    %   Laplacian and biharmonic values those of a quartic polynomial. The
    %   values are checked at an interior point of every polygon. The projection onto the discontinuous polynomials
    %   on the subtriangulation of the polygons reproduces polynomials, also
    %   at the tip of a slit and for a given discontinuous element. For an
    %   elementwise scaled linear function, ErrorForm gives the norms of the
    %   jumps of the function and of its normal derivative over the interior
    %   edges, and the squared EstimatorForm of every polygon is half the
    %   sum over its interior edges.
    %
    %   See also PDGk, PDGkVector, PDGkMatrix, PDGkSymMatrix, FESymMatrix,
    %   MeshPoly, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        polyMesh = struct( ...
            'SquareT',struct('class',@MeshPolyT,'name','Square','type','T'), ...
            'LshapeQ',struct('class',@MeshPolyQ,'name','Lshape','type','Q'), ...
            'SquareH',struct('class',@MeshPolyH,'name','Square','type','H'))
        subMesh = struct( ...
            'SquareH',struct('class',@MeshPolyH,'name','Square','type','H'), ...
            'LshapeH',struct('class',@MeshPolyH,'name','Lshape','type','H'), ...
            'SlitH',struct('class',@MeshPolyH,'name','Slit','type','H'), ...
            'SlitQ',struct('class',@MeshPolyQ,'name','Slit','type','Q'))
        element = struct( ...
            'PDG3',struct('space',@(m)PDGk(3,m),'field','scalar2d', ...
            'value',@(F)F.v,'grad',@(F)F.g), ...
            'PDG6',struct('space',@(m)PDGk(6,m),'field','sextic2d', ...
            'value',@(F)F.v,'grad',@(F)F.g), ...
            'PDGVector2',struct('space',@(m)PDGkVector(2,m),'field','vector2d', ...
            'value',@(F)F.v,'grad',@(F)FormBasisTest.flat(F.J)), ...
            'PDGMatrix1',struct('space',@(m)PDGkMatrix(1,m),'field','matrix2d', ...
            'value',@(F)F.s,'grad',[]), ...
            'PDGSymMatrix1',struct('space',@(m)PDGkSymMatrix(1,m),'field','symmatrix2d', ...
            'value',@(F)F.s,'grad',[]))
        derivative = struct( ...
            'D2_PDG6',struct('method','getD2Basis','L',@(F)F.D2), ...
            'D3_PDG6',struct('method','getD3Basis','L',@(F)F.D3,'relTol',1e-8), ...
            'D4_PDG6',struct('method','getD4Basis','L',@(F)F.D4,'relTol',1e-7))
    end

    methods (Test)

        function testReproducesPolynomial(testCase,polyMesh,element)
            mesh = PolygonalElementTest.create(polyMesh);
            fe = element.space(mesh);
            pts = PolygonalElementTest.interiorPoints(mesh);

            x = FormValueTest.coefficients(fe,element.field);

            F = FormValueTest.evalField(element.field,pts);
            testCase.verifyEqual(PolygonalElementTest.contract(fe,'getBasis',pts,x), ...
                element.value(F),'AbsTol',testCase.tol,'values');
            if ~isempty(element.grad)
                testCase.verifyEqual(PolygonalElementTest.contract(fe,'getGradBasis',pts,x), ...
                    element.grad(F),'AbsTol',testCase.tol,'gradients');
            end
        end

        function testHigherDerivatives(testCase,polyMesh,derivative)
            mesh = PolygonalElementTest.create(polyMesh);
            fe = PDGk(6,mesh);
            pts = PolygonalElementTest.interiorPoints(mesh);
            x = FormValueTest.coefficients(fe,'sextic2d');

            val = PolygonalElementTest.contract(fe,derivative.method,pts,x);

            expected = derivative.L(FormValueTest.evalField('sextic2d',pts));
            % the higher derivatives on the small polygons lose digits
            testCase.verifyEqual(val,expected, ...
                'AbsTol',DerivativeBasisTest.relTol(testCase,derivative)*max(1,max(abs(expected(:)))), ...
                'derivatives of the polynomial');
        end

        function testHigherDerivativeValues(testCase,polyMesh)
            mesh = PolygonalElementTest.create(polyMesh);
            fe = PDGk(4,mesh);
            pts = PolygonalElementTest.interiorPoints(mesh);
            parts = 1:mesh.nrElems;
            x = FormValueTest.coefficients(fe,'quartic2d');

            F = FormValueTest.evalField('quartic2d',pts);
            testCase.verifyEqual(fe.getD2Value(pts,[],parts,x),F.D2, ...
                'AbsTol',testCase.tol,'Hessian');
            testCase.verifyEqual(fe.getLaplaceValue(pts,[],parts,x),sum(F.D2(:,1:2),2), ...
                'AbsTol',testCase.tol,'Laplacian');
            testCase.verifyEqual(fe.getBiharmonicValue(pts,[],parts,x),F.bih, ...
                'AbsTol',1e-7*max(abs(F.bih)),'biharmonic');
        end

        function testElementwisePolynomial(testCase,polyMesh,element)
            mesh = PolygonalElementTest.create(polyMesh);
            fe = element.space(mesh);
            pts = PolygonalElementTest.interiorPoints(mesh);
            scale = 1 + (1:mesh.nrElems)'/mesh.nrElems;

            x = fe.project(@(z,zref,ind)scale(ind(:)).*FormValueTest.fieldValue(element.field,z,'s'));

            F = FormValueTest.evalField(element.field,pts);
            testCase.verifyEqual(PolygonalElementTest.contract(fe,'getBasis',pts,x), ...
                scale.*element.value(F),'AbsTol',testCase.tol);
        end

        function testProjectSubmeshDG(testCase,subMesh)
            mesh = subMesh.class(Geometry(struct('name',subMesh.name,'type',subMesh.type)));
            fe = PDGk(2,mesh);
            f = @(z)(1 + z(:,1) - 2*z(:,2)).^2;
            x = fe.project(f);

            [c,dgfe] = fe.projectSubmeshDG(x);
            c2 = fe.projectSubmeshDG(x,dgfe);

            ref = [0.2 0.3];
            val = dgfe.getValue([],ref,[],c);
            testCase.verifyGreaterThan(dgfe.mesh.area4e,0,'positive areas');
            testCase.verifyEqual(val(:),f(dgfe.mesh.point4e(ref)),'AbsTol',testCase.tol, ...
                'reproduces the polynomial');
            testCase.verifyEqual(c2,c,'AbsTol',testCase.tol,'given element');
        end

        function testJumpForms(testCase,subMesh)
            mesh = subMesh.class(Geometry(struct('name',subMesh.name,'type',subMesh.type)));
            fe = PDGk(1,mesh);
            s = (1:mesh.nrElems)'/mesh.nrElems;
            x = fe.project(@(z,zref,ind)(1 + s(ind(:))).*(z*[1;2]));

            err = ErrorForm('([u])^2 dsi').assemble(x,fe);
            errGrad = ErrorForm('([grad(u)]*n)^2 dsi').assemble(x,fe);
            eta4e = EstimatorForm('([u])^2 dsi').assemble(x,fe);

            % the jump of (1 + s_T)(x + 2y) is (s_1 - s_2)(x + 2y), integrated
            % exactly by the 3-point Gauss rule
            e4ed = mesh.e4ed;
            inner = e4ed(:,2) > 0;
            jump = zeros(mesh.nrEdges,1);
            jump(inner) = s(e4ed(inner,1)) - s(e4ed(inner,2));
            length4ed = mesh.length4ed;
            P1 = mesh.c4n(mesh.n4ed(:,1),:);
            P2 = mesh.c4n(mesh.n4ed(:,2),:);
            gauss = 0.5 + [-1 0 1]*sqrt(15)/10;
            weights = [5 8 5]/18;
            side = zeros(mesh.nrEdges,1);
            for q = 1 : 3
                side = side + weights(q)*length4ed.*(jump.*((P1 + gauss(q)*(P2 - P1))*[1;2])).^2;
            end
            normal4ed = mesh.normal4ed;
            sideGrad = length4ed.*(jump.*(normal4ed*[1;2])).^2;
            testCase.verifyEqual(err,sqrt(sum(side(inner))),'AbsTol',testCase.tol,'jump');
            testCase.verifyEqual(errGrad,sqrt(sum(sideGrad(inner))),'AbsTol',testCase.tol, ...
                'jump of the normal derivative');
            testCase.verifyEqual(eta4e.^2,accumarray([e4ed(inner,1);e4ed(inner,2)], ...
                [side(inner);side(inner)]/2,[mesh.nrElems,1]),'AbsTol',testCase.tol,'estimator');
        end

        function testSymMatrixOperators(testCase,polyMesh)
            mesh = PolygonalElementTest.create(polyMesh);
            fe = PDGkSymMatrix(1,mesh);
            pts = PolygonalElementTest.interiorPoints(mesh);
            parts = 1:mesh.nrElems;
            x = FormValueTest.coefficients(fe,'symmatrix2d');

            div = fe.getDivValue(pts,[],parts,x);
            rot = fe.getRotValue(pts,[],parts,x);
            tr = fe.getTrValue(pts,[],parts,x);
            divBasis = PolygonalElementTest.contract(fe,'getDivBasis',pts,x);

            F = FormValueTest.evalField('symmatrix2d',pts);
            testCase.verifyEqual(div,F.div,'AbsTol',testCase.tol,'row divergence');
            testCase.verifyEqual(divBasis,F.div,'AbsTol',testCase.tol,'divergence basis');
            testCase.verifyEqual(rot,F.rowRot,'AbsTol',testCase.tol,'row-wise rot');
            testCase.verifyEqual(tr,FormBasisTest.tr(F.M),'AbsTol',testCase.tol,'trace');
        end

    end

    methods (Static)

        function mesh = create(polyMesh)
            %create - Polygonal mesh, refined once at the first polygon
            mesh = polyMesh.class(Geometry(struct('name',polyMesh.name,'type',polyMesh.type)));
            r4e = false(mesh.nrElems,1);
            r4e(1) = true;
            mesh = mesh.refine(r4e);
        end

        function pts = interiorPoints(mesh)
            %interiorPoints - Points between the centroid and a vertex
            first = cellfun(@(nodes)nodes(1),mesh.n4e);
            pts = 0.7*mesh.mid4e + 0.3*mesh.c4n(first,:);
        end

        function val = contract(fe,method,pts,x)
            %contract - Contract a basis at PTS with the coefficients X
            nrElems = fe.mesh.nrElems;
            basis = fe.(method)(pts,[],1:nrElems);
            val = reshape(sum(basis.*x(fe.dof4e),2),nrElems,[]);
        end

    end

end
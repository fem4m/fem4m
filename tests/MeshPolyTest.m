classdef MeshPolyTest < UnitTestCase
    %MeshPolyTest - Unit tests of the polygonal meshes
    %   The triangle (T), quadrilateral (Q) and hexagonal (H) geometries of
    %   the square, the L-shaped domain and the slit domain give polygonal
    %   meshes whose areas add up to the area of the domain, with the
    %   diameter, centroid, edge lengths, tangents and normals of the
    %   polygons, and whose boundary has the length of the boundary of the
    %   triangulation of the domain. The subtriangulation at the centroids,
    %   or at a vertex for polygons not star-shaped to their centroid such
    %   as at the tip of a slit, is a conforming triangulation with
    %   consistent parent elements and edges. Refinement of MeshPolyT, MeshPolyQ and MeshPolyH preserves
    %   these properties, graded refinement meets its grading criterion,
    %   and the union of two meshes merges the common nodes and makes the
    %   common side interior, whatever its boundary conditions in the two
    %   meshes. Refinement of a circle projects the boundary nodes onto the
    %   circle. The symmetry boundary n4Sb keeps its fixed axis and,
    %   together with the Dirichlet and Neumann boundary, covers the
    %   boundary of the Dirichlet mesh refined in the same way; it survives
    %   the union of meshes, saving and loading, and the removal of ghost
    %   nodes. The boundary lists of the removal of ghost nodes are
    %   optional. A slit of height slitdelta keeps the fixed axes of the
    %   slit of height 0. Saved meshes reload with all digits, also in
    %   the subclasses (with hanging nodes, and with the files of the
    %   type first), meshes from node lists alone get a Dirichlet
    %   boundary, and one grading exponent serves several points.
    %
    %   See also MeshPoly, MeshPolyT, MeshPolyQ, MeshPolyH, Geometry,
    %   UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        polyMesh = struct( ...
            'SquareT',struct('class',@MeshPolyT,'name','Square','type','T','area',1), ...
            'SquareQ',struct('class',@MeshPolyQ,'name','Square','type','Q','area',1), ...
            'SquareH',struct('class',@MeshPolyH,'name','Square','type','H','area',1), ...
            'SquareHPoly',struct('class',@MeshPoly,'name','Square','type','H','area',1), ...
            'LshapeT',struct('class',@MeshPolyT,'name','Lshape','type','T','area',3), ...
            'LshapeQ',struct('class',@MeshPolyQ,'name','Lshape','type','Q','area',3), ...
            'LshapeH',struct('class',@MeshPolyH,'name','Lshape','type','H','area',3), ...
            'SlitT',struct('class',@MeshPolyT,'name','Slit','type','T','area',4), ...
            'SlitQ',struct('class',@MeshPolyQ,'name','Slit','type','Q','area',4), ...
            'SlitH',struct('class',@MeshPolyH,'name','Slit','type','H','area',4))
        refinable = struct( ...
            'SquareT',struct('class',@MeshPolyT,'name','Square','type','T','area',1), ...
            'SquareQ',struct('class',@MeshPolyQ,'name','Square','type','Q','area',1), ...
            'SquareH',struct('class',@MeshPolyH,'name','Square','type','H','area',1), ...
            'LshapeT',struct('class',@MeshPolyT,'name','Lshape','type','T','area',3), ...
            'LshapeQ',struct('class',@MeshPolyQ,'name','Lshape','type','Q','area',3), ...
            'LshapeH',struct('class',@MeshPolyH,'name','Lshape','type','H','area',3))
        % aligned: the geometry numbers the fixed axis (1 or 2) instead of 0
        symmetric = struct( ...
            'SquareT',struct('class',@MeshPolyT,'name','Square','type','T','bdry','DSNS','aligned',true,'area',1), ...
            'SquareQ',struct('class',@MeshPolyQ,'name','Square','type','Q','bdry','SSSS','aligned',true,'area',1), ...
            'SquareH',struct('class',@MeshPolyH,'name','Square','type','H','bdry','SDSN','aligned',true,'area',1), ...
            'SquareHAll',struct('class',@MeshPolyH,'name','Square','type','H','bdry','SSSS','aligned',true,'area',1), ...
            'LshapeT',struct('class',@MeshPolyT,'name','Lshape','type','T','bdry','NS','aligned',true,'area',3), ...
            'LshapeQ',struct('class',@MeshPolyQ,'name','Lshape','type','Q','bdry','SN','aligned',true,'area',3), ...
            'LshapeH',struct('class',@MeshPolyH,'name','Lshape','type','H','bdry','SN','aligned',true,'area',3), ...
            'SlitQ',struct('class',@MeshPolyQ,'name','Slit','type','Q','bdry','SNS','aligned',true,'area',4), ...
            'SlitH',struct('class',@MeshPolyH,'name','Slit','type','H','bdry','SNS','aligned',true,'area',4))
        % boundary conditions of the unit square and the unit square right
        % of it, whose common side is the side 2 of the left square and
        % the side 4 of the right one
        plusBdry = struct( ...
            'DD',struct('left','DDDD','right','DDDD'), ...
            'NDDD',struct('left','NDDD','right','NDDD'), ...
            'SS',struct('left','DSDD','right','SDSS'), ...
            'DN',struct('left','DDDD','right','NNNN'), ...
            'SD',struct('left','SSNN','right','DDDD'), ...
            'NS',struct('left','DNDS','right','NSNS'))
        circle = struct('D','D','DSDS','DSDS','SSSS','SSSS','NSNS','NSNS')
    end

    methods (Test)

        function testGeometry(testCase,polyMesh)
            mesh = MeshPolyTest.create(polyMesh);

            MeshPolyTest.verifyPolyMesh(testCase,mesh,polyMesh);
        end

        function testSubtriangulateMid(testCase,polyMesh)
            mesh = MeshPolyTest.create(polyMesh);

            MeshPolyTest.verifySubtriangulation(testCase,mesh,polyMesh.area);
        end

        function testRefine(testCase,refinable)
            mesh = MeshPolyTest.create(refinable);
            nrElems = mesh.nrElems;
            r4e = false(nrElems,1);
            r4e(1:2:end) = true;

            refined = mesh.refine(r4e);

            testCase.verifyGreaterThan(refined.nrElems,nrElems);
            MeshPolyTest.verifyPolyMesh(testCase,refined,refinable);
            MeshPolyTest.verifySubtriangulation(testCase,refined,refinable.area);
            % a second refinement of the first element
            r4e = false(refined.nrElems,1);
            r4e(1) = true;
            twice = refined.refine(r4e);
            MeshPolyTest.verifyPolyMesh(testCase,twice,refinable);
            MeshPolyTest.verifySubtriangulation(testCase,twice,refinable.area);
        end

        function testRefineKeepsTriangles(testCase)
            % MeshPolyH keeps marked triangles, the boundary edges of the
            % unmarked triangles must not be listed twice
            c4n = [0 0; 1 0; 1 1; 0 1; 2 0; 0.5 2; -1 0];
            n4e = {[1 2 3 4]; [2 5 3]; [3 5 6]; [4 3 6]; [7 1 4]; [7 4 6]};
            mesh = MeshPolyH(c4n,n4e,[7 1; 1 2; 2 5],[5 6; 6 7]);

            for r4e = {[true;true;false;false;false;false],[0;1;0;0;0;0]}
                refined = mesh.refine(r4e{1});

                length4ed = refined.length4ed;
                testCase.verifyEqual(sort([refined.ed4Db;refined.ed4Nb]), ...
                    find(refined.e4ed(:,2) == 0),'boundary edges');
                testCase.verifyEqual(sum(length4ed(refined.ed4Db)),3,'AbsTol',testCase.tol);
                testCase.verifyEqual(sum(length4ed(refined.ed4Nb)),5,'AbsTol',testCase.tol);
                MeshPolyTest.verifySubtriangulation(testCase,refined,3);
            end
        end

        function testClosure(testCase)
            mesh = MeshPolyTest.create(struct('class',@MeshPolyH,'name','Lshape','type','H'));
            r4e = false(mesh.nrElems,1);
            r4e(1) = true;

            closed = mesh.closure(r4e);

            testCase.verifyTrue(all(closed(r4e)),'marked elements stay marked');
            testCase.verifyEqual(mesh.closure(closed),closed,'idempotent');
        end

        function testRefineGraded(testCase,refinable)
            mesh = MeshPolyTest.create(refinable);
            mu = 0.5;
            H = 0.5;

            refined = mesh.refineGraded(mu,H,[0 0]);

            hT = refined.diam4e;
            rT = sqrt(sum(refined.mid4e.^2,2));
            rT(rT < eps) = hT(rT < eps)/2;
            testCase.verifyLessThanOrEqual(hT,H + testCase.tol,'mesh size');
            testCase.verifyLessThanOrEqual(hT,H*rT.^(1-mu) + testCase.tol,'grading');
            MeshPolyTest.verifyPolyMesh(testCase,refined,refinable);
        end

        function testPlus(testCase)
            left = MeshPoly(Geometry(struct('name','Square','type','Q')));
            right = MeshPoly(Geometry(struct('name','Square','type','Q','midpoint',[1.5 0.5])));
            shared = nnz(ismember(left.c4n,right.c4n,'rows'));

            union = left + right;

            testCase.verifyEqual(union.nrElems,left.nrElems + right.nrElems);
            testCase.verifyEqual(union.nrNodes,left.nrNodes + right.nrNodes - shared);
            testCase.verifyEqual(sum(union.area4e),2,'AbsTol',testCase.tol);
            length4ed = union.length4ed;
            testCase.verifyEqual(sum(length4ed(union.ed4Db)),6,'AbsTol',testCase.tol, ...
                'the common side is interior');
        end

        function testPlusBoundary(testCase,plusBdry)
            left = MeshPoly(Geometry(struct('name','Square','type','Q','bdry',plusBdry.left)));
            right = MeshPoly(Geometry(struct('name','Square','type','Q','bdry',plusBdry.right,'midpoint',[1.5 0.5])));

            union = left + right;

            testCase.verifyEqual(union.nrEdges,left.nrEdges + right.nrEdges - 1);
            testCase.verifyEqual(sort([union.ed4Db;union.ed4Nb;union.ed4Sb]),find(union.e4ed(:,2) == 0), ...
                'the common side is interior');
            % the other sides keep their boundary conditions
            kept = [plusBdry.left([1 3 4]),plusBdry.right(1:3)];
            length4ed = union.length4ed;
            testCase.verifyEqual(sum(length4ed(union.ed4Db)),nnz(kept == 'D'),'AbsTol',testCase.tol);
            testCase.verifyEqual(sum(length4ed(union.ed4Nb)),nnz(kept == 'N'),'AbsTol',testCase.tol);
            testCase.verifyEqual(sum(length4ed(union.ed4Sb)),nnz(kept == 'S'),'AbsTol',testCase.tol);
            if ~isempty(union.n4Sb)
                tangents = union.tangent4ed;
                testCase.verifyEqual(union.n4Sb(:,3),1 + double(abs(tangents(union.ed4Sb,2)) < testCase.tol), ...
                    'fixed axes');
            end
        end

        function testRemoveGhostNodes(testCase)
            % node 1 is a ghost node, the boundary lists are optional
            c4n = [9 9; 0 0; 1 0; 0 1];
            n4e = {[2 3 4]};

            [c4nRed,n4eRed] = MeshPoly.removeGhostNodes(c4n,n4e);
            [~,~,n4Db] = MeshPoly.removeGhostNodes(c4n,n4e,[3 4]);

            testCase.verifyEqual(c4nRed,[0 0; 1 0; 0 1]);
            testCase.verifyEqual(n4eRed,{[1 2 3]});
            testCase.verifyEqual(n4Db,[2 3]);
        end

        function testColor4nSingleElement(testCase)
            % the nodes of one element get different colors on every level
            mesh = MeshPolyQ(Geometry(struct('name','Square','type','Q')));
            for levels = 1 : 2
                testCase.verifyEqual(sort(mesh.computeColor4n(levels)),(1:4)', ...
                    sprintf('levels %d',levels));
            end
        end

        function testSymmetryBoundary(testCase,symmetric)
            mesh = MeshPolyTest.create(symmetric);
            dirichlet = MeshPolyTest.create(rmfield(symmetric,'bdry'));

            MeshPolyTest.verifySymmetryBoundary(testCase,mesh,dirichlet,symmetric);
            MeshPolyTest.verifySubtriangulation(testCase,mesh,symmetric.area);
            for k = 1 : 2
                r4e = false(mesh.nrElems,1);
                r4e(k:2:end) = true;
                mesh = mesh.refine(r4e);
                dirichlet = dirichlet.refine(r4e);
                MeshPolyTest.verifyPolyMesh(testCase,mesh,symmetric);
                MeshPolyTest.verifySymmetryBoundary(testCase,mesh,dirichlet,symmetric);
                MeshPolyTest.verifySubtriangulation(testCase,mesh,symmetric.area);
            end
        end

        function testSymmetryPlus(testCase)
            left = MeshPoly(Geometry(struct('name','Square','type','Q','bdry','S')));
            right = MeshPoly(Geometry(struct('name','Square','type','Q','bdry','S','midpoint',[1.5 0.5])));

            union = left + right;

            length4ed = union.length4ed;
            testCase.verifyEqual(sum(length4ed(union.ed4Sb)),6,'AbsTol',testCase.tol, ...
                'the common side is interior');
            testCase.verifyEqual(sort(union.ed4Sb),find(union.e4ed(:,2) == 0));
            testCase.verifyEqual(union.n4ed(union.ed4Sb,:),union.n4Sb(:,1:2));
            tangents = union.tangent4ed;
            testCase.verifyEqual(union.n4Sb(:,3),1 + double(abs(tangents(union.ed4Sb,2)) < testCase.tol), ...
                'fixed axes');
        end

        function testSymmetrySlitDelta(testCase)
            % the slit of height slitdelta keeps the fixed axes of the slit
            % of height 0, also on the edges at its tip that are not
            % horizontal
            for type = 'TQH'
                slit = Geometry(struct('name','Slit','type',type,'bdry','SNS'));
                wide = Geometry(struct('name','Slit','type',type,'bdry','SNS','slitdelta',0.1));

                testCase.verifyEqual(wide.n4Sb,slit.n4Sb,type);
                testCase.verifyTrue(all(wide.n4Sb(:,3) > 0),type);
            end
        end

        function testSymmetrySaveLoad(testCase)
            folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            mesh = MeshPolyH(Geometry(struct('name','Square','type','H','bdry','SNSD')));
            prefix = fullfile(folder,'mesh');

            mesh.save(prefix);
            loaded = MeshPoly(prefix);
            subclass = MeshPolyH(prefix);

            testCase.verifyEqual(loaded.c4n,mesh.c4n,'all digits');
            testCase.verifyEqual(loaded.n4e,mesh.n4e);
            testCase.verifyEqual(subclass.n4e,mesh.n4e,'MeshPolyH reads the saved files');
            testCase.verifyEqual(loaded.n4Sb,mesh.n4Sb);
            testCase.verifyEqual(loaded.ed4Sb,mesh.ed4Sb);
            testCase.verifyEqual(loaded.SbEd4n,mesh.SbEd4n);
        end

        function testHangingNodesSaveLoad(testCase)
            % the underlying quadrilaterals and triangles and the hanging
            % nodes survive save and load
            folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            meshes = {MeshPolyQ(Geometry(struct('name','Square','type','Q','bdry','DSNS'))), ...
                MeshPolyT(Geometry(struct('name','Square','type','T','bdry','DSNS')))};
            nrVertices = [4 3];
            for k = 1 : numel(meshes)
                mesh = meshes{k};
                for level = 1 : 3
                    r4e = false(mesh.nrElems,1);
                    r4e(1) = true;
                    mesh = mesh.refine(r4e);
                end
                testCase.assertGreaterThan(max(cellfun(@numel,mesh.n4e)),nrVertices(k),'hanging nodes');
                prefix = fullfile(folder,sprintf('mesh%d',k));

                mesh.save(prefix);
                loaded = feval(class(mesh),prefix);

                testCase.verifyEqual(loaded.n4e,mesh.n4e,class(mesh));
                testCase.verifyEqual(loaded.hnodes,mesh.hnodes,class(mesh));
                r4e = true(mesh.nrElems,1);
                testCase.verifyEqual(loaded.refine(r4e).n4e,mesh.refine(r4e).n4e,class(mesh));
            end
        end

        function testTypedFilesFirst(testCase)
            % MeshPolyQ(NAME) reads NAMEQ_*.dat before NAME_*.dat
            folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            prefix = fullfile(folder,'base');
            MeshPoly(Geometry(struct('name','Lshape','type','T'))).save(prefix);
            quad = MeshPolyQ(Geometry(struct('name','Square','type','Q')));
            quad.save([prefix,'Q']);

            loaded = MeshPolyQ(prefix);

            testCase.verifyEqual(loaded.n4quad,quad.n4quad);
        end

        function testGradedScalarExponent(testCase,refinable)
            % one exponent for several grading points
            points = [0 0;0.5 0.5];
            scalar = MeshPolyTest.create(refinable).refineGraded(0.5,0.3,points);
            vector = MeshPolyTest.create(refinable).refineGraded([0.5 0.5],0.3,points,[1 1]);

            testCase.verifyEqual(scalar.n4e,vector.n4e);
            testCase.verifyEqual(scalar.c4n,vector.c4n);
        end

        function testNumericConstructor(testCase)
            geom = Geometry(struct('name','Square','type','Q'));

            % all boundary edges are Dirichlet edges
            mesh = MeshPoly(geom.c4n,geom.n4e);
            quad = MeshPolyQ(geom.c4n,geom.n4e);

            testCase.verifyEqual(sortrows(sort(mesh.n4Db,2)),sortrows(sort(geom.n4Db,2)));
            testCase.verifyEmpty(mesh.n4Nb);
            testCase.verifyEqual(quad.n4Db,mesh.n4Db,'MeshPolyQ');
        end

        function testSymmetryRemoveGhostNodes(testCase)
            % node 1 is a ghost node, the axis 2 of the symmetry boundary
            % must not be renumbered as a node
            c4n = [9 9; 0 0; 1 0; 1 1; 0 1];
            n4e = {[2 3 4 5]};

            [c4n,n4e,n4Db,n4Nb,n4Sb] = MeshPoly.removeGhostNodes(c4n,n4e, ...
                [5 2],[3 4],[2 3 2; 4 5 2]);

            testCase.verifyEqual(c4n,[0 0; 1 0; 1 1; 0 1]);
            testCase.verifyEqual(n4e,{[1 2 3 4]});
            testCase.verifyEqual(n4Db,[4 1]);
            testCase.verifyEqual(n4Nb,[2 3]);
            testCase.verifyEqual(n4Sb,[1 2 2; 3 4 2]);
        end

        function testSymmetryDisplay(testCase)
            mesh = MeshPolyQ(Geometry(struct('name','Square','type','Q','bdry','DSNS')));
            fig = figure('Visible','off');
            testCase.addTeardown(@close,fig);

            text = evalc('disp(mesh)');
            mesh.plotBd;

            testCase.verifySubstring(text,'nodes for symmetric boundary');
            symmetry = findobj(fig,'Type','patch','LineStyle',':');
            testCase.verifyNumElements(symmetry,1);
            testCase.verifySize(symmetry.XData,[2 size(mesh.n4Sb,1)]);
        end

        function testCircle(testCase,circle)
            mesh = Mesh(Geometry(struct('name','Circle','bdry',circle)));

            refined = mesh.refine('uniform',2);

            n4bd = refined.n4ed(refined.e4ed(:,2) == 0,:);
            r = sqrt(sum(refined.c4n(unique(n4bd),:).^2,2));
            testCase.verifyEqual(r,ones(size(r)),'AbsTol',testCase.tol);
            testCase.verifyEqual(sum(refined.area4e),pi,'RelTol',0.05);
        end

    end

    methods (Static)

        function mesh = create(polyMesh)
            %create - Polygonal mesh of the class and geometry of POLYMESH
            opts = struct('name',polyMesh.name,'type',polyMesh.type);
            if isfield(polyMesh,'bdry')
                opts.bdry = polyMesh.bdry;
            end
            mesh = polyMesh.class(Geometry(opts));
        end

        function verifyPolyMesh(testCase,mesh,polyMesh)
            %verifyPolyMesh - Verify areas, diameters, centroids and edges
            tol = testCase.tol;
            area4e = mesh.area4e;
            diam4e = mesh.diam4e;
            mid4e = mesh.mid4e;
            tangent4e = mesh.tangent4e;
            normal4e = mesh.normal4e;
            testCase.verifyGreaterThan(area4e,0,'positive area');
            testCase.verifyEqual(sum(area4e),polyMesh.area,'AbsTol',tol,'total area');
            testCase.verifyEqual(mesh.nrN4e,cellfun(@numel,mesh.n4e(:)));
            for T = 1 : mesh.nrElems
                P = mesh.c4n(mesh.n4e{T},:);
                shape = polyshape(P,'Simplify',false);
                testCase.verifyEqual(area4e(T),area(shape),'AbsTol',tol,sprintf('area4e(%d)',T));
                testCase.verifyEqual(diam4e(T),max(vecnorm(P - permute(P,[3 2 1]),2,2),[],'all'),'AbsTol',tol,sprintf('diam4e(%d)',T));
                if isprop(mesh,'n4tri')
                    expected = mean(mesh.c4n(mesh.n4tri(T,:),:),1);
                elseif isprop(mesh,'n4quad')
                    expected = mean(mesh.c4n(mesh.n4quad(T,:),:),1);
                else
                    [x,y] = centroid(shape);
                    expected = [x,y];
                end
                testCase.verifyEqual(mid4e(T,:),expected,'AbsTol',tol,sprintf('mid4e(%d)',T));
                % tangents and outer normals of the polygon
                t = reshape(tangent4e{T},[],2);
                n = reshape(normal4e{T},[],2);
                edges = P([2:end,1],:) - P;
                lengths = sqrt(sum(edges.^2,2));
                testCase.verifyEqual(t,edges./lengths,'AbsTol',tol,sprintf('tangent4e{%d}',T));
                testCase.verifyEqual(n,[t(:,2),-t(:,1)],'AbsTol',tol,sprintf('normal4e{%d}',T));
                testCase.verifyEqual(sum(lengths.*n,1),[0 0],'AbsTol',tol, ...
                    sprintf('closed polygon %d',T));
                % counterclockwise vertices, so that (t2,-t1) points outward
                % also for the nonconvex elements at reentrant corners
                Q = P([2:end,1],:);
                testCase.verifyGreaterThan(sum(P(:,1).*Q(:,2) - Q(:,1).*P(:,2)),0, ...
                    sprintf('counterclockwise element %d',T));
            end

            % edges
            length4ed = mesh.length4ed;
            tangent4ed = mesh.tangent4ed;
            normal4ed = mesh.normal4ed;
            edges = mesh.c4n(mesh.n4ed(:,2),:) - mesh.c4n(mesh.n4ed(:,1),:);
            testCase.verifyEqual(length4ed,sqrt(sum(edges.^2,2)),'AbsTol',tol);
            testCase.verifyEqual(tangent4ed.*length4ed,edges,'AbsTol',tol);
            testCase.verifyEqual(sum(normal4ed.*tangent4ed,2),zeros(mesh.nrEdges,1),'AbsTol',tol);
            testCase.verifyEqual(sum(normal4ed.^2,2),ones(mesh.nrEdges,1),'AbsTol',tol);
            testCase.verifyEqual(mesh.point4ed(1/2),(mesh.c4n(mesh.n4ed(:,1),:) ...
                + mesh.c4n(mesh.n4ed(:,2),:))/2,'AbsTol',tol);
            edNr4ed = mesh.edNr4ed;
            for i = 1 : mesh.nrEdges
                for side = 1 : 2
                    T = mesh.e4ed(i,side);
                    if T > 0
                        testCase.verifyEqual(mesh.ed4e{T}(edNr4ed(i,side)),i, ...
                            sprintf('edNr4ed(%d,%d)',i,side));
                    end
                end
            end

            % the boundary has the length of the boundary of the domain
            triangulation = Mesh(Geometry(struct('name',polyMesh.name)));
            triLength4ed = triangulation.length4ed;
            boundary = mesh.e4ed(:,2) == 0;
            testCase.verifyEqual(sort([mesh.ed4Db;mesh.ed4Nb;mesh.ed4Sb]),find(boundary),'boundary edges');
            testCase.verifyEqual(sum(length4ed(boundary)), ...
                sum(triLength4ed(triangulation.e4ed(:,2) == 0)),'AbsTol',tol,'boundary length');
        end

        function verifySymmetryBoundary(testCase,mesh,dirichlet,polyMesh)
            %verifySymmetryBoundary - Verify the symmetry boundary data
            %   The mesh DIRICHLET of the same geometry with Dirichlet
            %   boundary only is refined in the same way, so its boundary
            %   is the union of all boundary parts of MESH.
            testCase.verifyEqual(mesh.c4n,dirichlet.c4n,'same refinement');
            testCase.verifyNotEmpty(mesh.n4Sb);
            testCase.verifyEqual(sortrows([mesh.n4Db;mesh.n4Nb;mesh.n4Sb(:,1:2)]), ...
                sortrows(dirichlet.n4Db),'boundary node pairs');
            testCase.verifyEqual(mesh.n4ed(mesh.ed4Sb,:),mesh.n4Sb(:,1:2),'n4ed of symmetry boundary');
            testCase.verifyEqual(mesh.SbEd4n(mesh.n4Sb(:,1),1),mesh.ed4Sb,'SbEd4n(:,1)');
            testCase.verifyEqual(mesh.SbEd4n(mesh.n4Sb(:,2),2),mesh.ed4Sb,'SbEd4n(:,2)');
            testCase.verifyEqual(nnz(mesh.SbEd4n),2*numel(mesh.ed4Sb),'SbEd4n of other nodes');
            % the y-coordinate is fixed (axis 2) on horizontal edges and the
            % x-coordinate (axis 1) on vertical ones
            if polyMesh.aligned
                tangents = mesh.tangent4ed;
                testCase.verifyEqual(mesh.n4Sb(:,3), ...
                    1 + double(abs(tangents(mesh.ed4Sb,2)) < testCase.tol),'fixed axes');
            else
                testCase.verifyEqual(mesh.n4Sb(:,3),zeros(size(mesh.n4Sb,1),1),'axes');
            end
        end

        function verifySubtriangulation(testCase,mesh,area)
            %verifySubtriangulation - Verify the triangulation at centroids
            area4e = mesh.area4e;
            length4ed = mesh.length4ed;

            [tri,parent4ed,children4ed] = mesh.subtriangulateMid;

            % one triangle per edge from the centroid, or two less from a
            % vertex for polygons not star-shaped to their centroid
            nrChildren = accumarray(tri.parent4e,1,[mesh.nrElems,1]);
            testCase.verifyTrue(all(nrChildren == mesh.nrN4e | nrChildren == mesh.nrN4e - 2), ...
                'triangles per polygon');
            MeshTest.verifyMesh(testCase,tri,tri,area);
            testCase.verifyEqual(accumarray(tri.parent4e,tri.area4e,[mesh.nrElems,1]), ...
                area4e,'AbsTol',testCase.tol,'areas of the children');
            triLength4ed = tri.length4ed;
            children = parent4ed > 0;
            testCase.verifyEqual(accumarray(parent4ed(children),triLength4ed(children), ...
                [mesh.nrEdges,1]),length4ed,'AbsTol',testCase.tol,'lengths of the child edges');
            testCase.verifyEqual(sort(children4ed),find(children));
            testCase.verifyTrue(issorted(parent4ed(children4ed)),'children4ed ordered by parent');
            testCase.verifyEqual(tri.n4Sb,mesh.n4Sb,'symmetry boundary');
        end

    end

end
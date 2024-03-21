classdef MeshTest < UnitTestCase
    %MeshTest - Unit tests of the refinement of Mesh and Mesh3d
    %   Uniform refinement of a triangular mesh splits every element into
    %   four, so that nrElems' = 4*nrElems, nrNodes' = nrNodes + nrEdges
    %   and nrEdges' = 2*nrEdges + 3*nrElems; uniform refinement of a
    %   tetrahedral mesh (three bisection sweeps) gives every element
    %   eight children, also for vertex colorings with more than four
    %   colors and for longest edge labels. One bisection sweep of the
    %   stored labels of the 3d geometries doubles the elements, and
    %   these labels are the longest edge labels. Longest edge labels
    %   give all four types of Arnold, Mukherjee, Pouly on a Delaunay
    %   mesh, refine the longest edge first, and refine single elements
    %   of every type, also with equally long edges, within their
    %   bisection tree. Coloring labels have type 0 and refine
    %   adaptively, also with permuted node ids. Inconsistent labels
    %   raise unmatchedFaces at the first refine that exposes them.
    %   Uniform and local (red-green-blue or newest vertex bisection)
    %   refinement and coarsening preserve the total area or volume, the
    %   length or area of the boundary and the Euler characteristic, keep
    %   all elements positively oriented, and yield conforming meshes with
    %   consistent edge, face and boundary data structures. Coarsening a
    %   uniform newest vertex bisection restores the initial boundary
    %   parts, including single boundary edges and the symmetry axes. Meshes
    %   survive save and load (with all digits of the coordinates, labels
    %   and symmetry axes, also for far-off rotated cubes) and large node
    %   ids. Zero uniform rounds give the identity refinement data, and
    %   graded refinement works on two elements and around several
    %   points.
    %   Invalid refinement options and refinement without marked edges
    %   issue warnings.
    %
    %   See also Mesh, Mesh3d, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        domain = struct( ...
            'Square',struct('name','Square','area',1), ...
            'SquareDNDN',struct('name','Square','bdry','DNDN','area',1), ...
            'Lshape',struct('name','Lshape','area',3), ...
            'Slit',struct('name','Slit','area',4), ...
            'DblConnected',struct('name','DblConnected','area',12))
        refinement = struct('RGB','RGB','NVB','NVB')
        bdryDomain = struct( ...
            'SquareDSNS',struct('name','Square','bdry','DSNS','area',1), ...
            'SquareSDNN',struct('name','Square','bdry','SDNN','area',1), ...
            'SquareRotated',struct('name','Square','bdry','SSSS','rotation',pi/6,'area',1), ...
            'LshapeSN',struct('name','Lshape','bdry','SN','area',3))
    end

    methods (Test)

        function testUniformRefinementCounts(testCase,domain)
            mesh = MeshTest.domainMesh(domain);

            refined = mesh.refine('uniform');

            testCase.verifyEqual(refined.nrElems,4*mesh.nrElems);
            testCase.verifyEqual(refined.nrNodes,mesh.nrNodes + mesh.nrEdges);
            testCase.verifyEqual(refined.nrEdges,2*mesh.nrEdges + 3*mesh.nrElems);
            testCase.verifyEqual(numel(refined.ed4Db),2*numel(mesh.ed4Db));
            testCase.verifyEqual(numel(refined.ed4Nb),2*numel(mesh.ed4Nb));
        end

        function testUniformRefinement(testCase,domain)
            mesh = MeshTest.domainMesh(domain);

            refined = mesh.refine('uniform',2);

            MeshTest.verifyMesh(testCase,refined,mesh,domain.area);
        end

        function testLocalRefinement(testCase,domain,refinement)
            mesh = MeshTest.domainMesh(domain).refine('uniform');
            marked4e = zeros(mesh.nrElems,1);
            marked4e(1:3:end) = 1;

            refined = mesh.refine(marked4e,refinement);

            testCase.verifyGreaterThan(refined.nrElems,mesh.nrElems);
            MeshTest.verifyMesh(testCase,refined,mesh,domain.area);
        end

        function testRepeatedLocalRefinement(testCase,domain,refinement)
            mesh = MeshTest.domainMesh(domain);

            refined = mesh;
            for k = 1 : 4
                marked4e = zeros(refined.nrElems,1);
                marked4e(1) = 1;
                refined = refined.refine(marked4e,refinement);
            end

            MeshTest.verifyMesh(testCase,refined,mesh,domain.area);
        end

        function testLogicalMarks(testCase,refinement)
            mesh = MeshTest.domainMesh(struct('name','Lshape')).refine('uniform');
            marked4e = zeros(mesh.nrElems,1);
            marked4e(1:3:end) = 1;
            marked4ed = zeros(mesh.nrEdges,1);
            marked4ed(1:5:end) = 1;

            % logical element and edge masks refine like 0/1 vectors
            testCase.verifyEqual(mesh.refine(logical(marked4e),refinement).n4e, ...
                mesh.refine(marked4e,refinement).n4e);
            testCase.verifyEqual(mesh.refine(logical(marked4ed),refinement).n4e, ...
                mesh.refine(marked4ed,refinement).n4e);
        end

        function testCoarsening(testCase,domain)
            mesh = MeshTest.domainMesh(domain).refine('uniform');
            marked4e = zeros(mesh.nrElems,1);
            marked4e(1) = 1;
            refined = mesh.refine(marked4e,'NVB');

            % every coarsening step undoes one bisection
            coarsened = refined.coarse(true(refined.nrElems,1));
            coarsened = coarsened.coarse(true(coarsened.nrElems,1));

            MeshTest.verifyMesh(testCase,coarsened,mesh,domain.area);
            testCase.verifyEqual(coarsened.nrElems,mesh.nrElems);
            testCase.verifyEqual(coarsened.nrNodes,mesh.nrNodes);
        end

        function testCoarseningBoundary(testCase,bdryDomain)
            mesh = MeshTest.domainMesh(bdryDomain);
            refined = mesh.refine('NVB','uniform');

            % two coarsening steps undo the uniform bisection, keeping
            % single boundary edges and the axes 0, 1, 2 of n4Sb
            coarsened = refined.coarse(true(refined.nrElems,1));
            coarsened = coarsened.coarse(true(coarsened.nrElems,1));

            MeshTest.verifyMesh(testCase,coarsened,mesh,bdryDomain.area);
            testCase.verifyEqual(coarsened.c4n,mesh.c4n);
            testCase.verifyEqual(sortrows(coarsened.n4Db),sortrows(mesh.n4Db),'n4Db');
            testCase.verifyEqual(sortrows(coarsened.n4Nb),sortrows(mesh.n4Nb),'n4Nb');
            testCase.verifyEqual(sortrows(coarsened.n4Sb),sortrows(mesh.n4Sb),'n4Sb');
        end

        function testInvalidRefinementOption(testCase)
            mesh = MeshTest.domainMesh(struct('name','Square'));

            refined = testCase.verifyWarning(@()mesh.refine('nonsense','uniform'), ...
                'fem4m:Mesh:unknownRefineOption');

            testCase.verifyEqual(refined.nrElems,4*mesh.nrElems);
        end

        function testNoMarkedEdges(testCase)
            mesh = MeshTest.domainMesh(struct('name','Square'));

            refined = testCase.verifyWarning(@()mesh.refine(zeros(mesh.nrElems,1)), ...
                'fem4m:Mesh:noMarkedEdges');

            testCase.verifyEqual(refined.nrElems,4*mesh.nrElems);
        end

        function testUniformRefinement3d(testCase)
            mesh = Mesh3d(Geometry(struct('name','Cube', ...
                'midpoint',[0.5,0.5,0.5],'length',1,'bdry','DNDNDN')));

            refined = mesh.refine('uniform');

            testCase.verifyEqual(refined.nrElems,8*mesh.nrElems);
            testCase.verifyEqual(refined.nrNodes,mesh.nrNodes + mesh.nrEdges);
            testCase.verifyEqual(numel(refined.f4Db),4*numel(mesh.f4Db));
            testCase.verifyEqual(numel(refined.f4Nb),4*numel(mesh.f4Nb));
            MeshTest.verifyMesh3d(testCase,refined,mesh,1);
            % N rounds give 8^N descendants per element
            uniform2 = mesh.refine('uniform',2);
            testCase.verifyEqual(accumarray(double(uniform2.parent4e(:)),1), ...
                64*ones(mesh.nrElems,1),'two rounds');
            testCase.verifyEqual(numel(uniform2.f4Db),16*numel(mesh.f4Db));
            MeshTest.verifyMesh3d(testCase,uniform2,mesh,1);
        end

        function testUniformColoring3d(testCase)
            % more than 4 colors: the three sweeps of a round need no
            % closure and every element gets 8 children per round
            c4n = [0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1; ...
                0.4 0.5 0.6;0.7 0.2 0.3;0.2 0.8 0.4;0.6 0.6 0.2;0.5 0.3 0.8];
            n4e = delaunay(c4n);
            mesh = Mesh3d(c4n,n4e,Mesh3d.getBoundary(n4e),[],'init','coloring');
            testCase.assertGreaterThan(max(mesh.color4n),4,'more than 4 colors');
            testCase.verifyTrue(all(mesh.t4e == 0 | mesh.t4e == 3),'type 0');

            refined = mesh.refine('uniform',2);

            testCase.verifyEqual(accumarray(double(refined.parent4e(:)),1), ...
                64*ones(mesh.nrElems,1));
            MeshTest.verifyMesh3d(testCase,refined,mesh,1);
        end

        function testZeroUniformRounds(testCase)
            mesh = MeshTest.domainMesh(struct('name','Square')).refine('uniform');

            unchanged = mesh.refine('uniform',0);

            testCase.verifyEqual(unchanged.n4e,mesh.n4e);
            testCase.verifyEqual(unchanged.parent4e,(1:mesh.nrElems)','identity parents');
            testCase.verifyEqual(full(unchanged.prolongation),eye(mesh.nrNodes),'identity P1');
            P = Pk(2,unchanged).prolongation(Pk(2,mesh));
            testCase.verifyEqual(full(P),eye(size(P,1)),'AbsTol',testCase.tol,'identity prolongation');
        end

        function testNumericConstructor(testCase)
            geom = Geometry(struct('name','Square'));

            % all boundary edges are Dirichlet edges
            mesh = Mesh(geom.c4n,geom.n4e);

            testCase.verifyEqual(sortrows(sort(mesh.n4Db,2)),sortrows(sort(geom.n4Db,2)));
            testCase.verifyEmpty(mesh.n4Nb);
        end

        function testSaveLoad(testCase)
            % empty lists delete the files of an earlier save
            mixed = Mesh(Geometry(struct('name','Square','bdry','DSNS')));
            dirichlet = Mesh(Geometry(struct('name','Square','bdry','DDDD')));
            file = [tempname,'_mesh'];
            cleanup = onCleanup(@()delete([file,'_*.dat']));

            mixed.save(file);
            loaded = Mesh(file);
            dirichlet.save(file);
            reloaded = Mesh(file);

            testCase.verifyEqual(loaded.n4Sb,mixed.n4Sb,'symmetry pairs with axis');
            testCase.verifyEqual(reloaded.n4Db,dirichlet.n4Db);
            testCase.verifyEmpty(reloaded.n4Nb,'no stale Neumann edges');
            testCase.verifyEmpty(reloaded.n4Sb,'no stale symmetry edges');
        end

        function testBisection3d(testCase)
            mesh = UnitTestCase.cubeMesh('DNDNDN');

            refined = mesh;
            for k = 1 : 3
                % refine marks by a 0/1 vector over the elements
                marked4e = zeros(refined.nrElems,1);
                marked4e([1,refined.nrElems]) = 1;
                refined = refined.refine(marked4e);
            end

            testCase.verifyGreaterThan(refined.nrElems,mesh.nrElems);
            MeshTest.verifyMesh3d(testCase,refined,mesh,1);
        end

        function testUniformAfterAdaptive3d(testCase)
            mesh = Mesh3d(Geometry(struct('name','Cube', ...
                'midpoint',[0.5,0.5,0.5],'length',1,'bdry','DNDNDN')));
            adapted = mesh.refine([true;false(mesh.nrElems-1,1)]);

            refined = adapted.refine('uniform');

            % the closure after the three sweeps completes the mesh
            testCase.verifyGreaterThanOrEqual(refined.nrElems,8*adapted.nrElems);
            MeshTest.verifyMesh3d(testCase,refined,mesh,1);
        end

        function testGradedRefinement(testCase)
            mesh = MeshTest.domainMesh(struct('name','Lshape'));
            mu = 0.5;
            H = 0.25;

            refined = mesh.refine('graded',mu,H,[0 0]);

            % no element violates the grading criterion at the corner
            length4ed = refined.length4ed;
            hT = max(length4ed(refined.ed4e),[],2);
            rT = sqrt(sum(refined.point4e([1 1]/3).^2,2));
            testCase.verifyLessThanOrEqual(hT,H + testCase.tol,'mesh size');
            testCase.verifyLessThanOrEqual(hT,H*rT.^(1-mu) + testCase.tol,'grading');
            testCase.verifyLessThan(min(hT),H/4,'refined towards the corner');
            MeshTest.verifyMesh(testCase,refined,mesh,3);
        end

        function testGradedBetaRefinement(testCase)
            mesh = MeshTest.domainMesh(struct('name','Lshape'));
            uniform = mesh.refine('uniform');
            mu = 2/3;
            R = 3/4;

            graded = mesh.refine('gradedbeta',mu,[0 0],R);
            plain = mesh.refine('gradedbeta',1,[0 0],R);

            testCase.verifyEqual(graded.n4e,uniform.n4e,'same topology');
            r = sqrt(sum(uniform.c4n.^2,2));
            moved = r > 0 & r < R;
            testCase.assertNotEmpty(find(moved,1));
            expected = uniform.c4n;
            expected(moved,:) = R*(r(moved)/R).^(1/mu).*uniform.c4n(moved,:)./r(moved);
            testCase.verifyEqual(graded.c4n,expected,'AbsTol',testCase.tol,'shifted nodes');
            testCase.verifyEqual(plain.c4n,uniform.c4n,'AbsTol',testCase.tol,'mu = 1');
            MeshTest.verifyMesh(testCase,graded,mesh,3);
        end

        function testGradedRedRefinement(testCase)
            mesh = MeshTest.domainMesh(struct('name','Lshape'));
            corner = find(all(mesh.c4n == 0,2));
            testCase.assertNumElements(corner,1);
            kappa = 0.25;

            red = mesh.refine('gradedred',0.5,[0 0]);
            graded = mesh.refine('gradedred',kappa,[0 0]);

            testCase.verifyEqual(red.nrElems,4*mesh.nrElems);
            testCase.verifyEqual(sortrows(red.c4n),sortrows(mesh.refine('uniform').c4n), ...
                'AbsTol',testCase.tol,'kappa = 1/2 is red refinement');
            % the new nodes on the edges at the corner divide them by kappa
            atCorner = find(any(mesh.n4ed == corner,2));
            other = sum(mesh.n4ed(atCorner,:),2) - corner;
            expected = kappa*mesh.c4n(other,:);
            newNode4ed = mesh.nrNodes + atCorner;
            testCase.verifyEqual(graded.c4n(newNode4ed,:),expected,'AbsTol',testCase.tol);
            testCase.verifyEqual(graded.parent4e,repmat((1:mesh.nrElems)',4,1));
            MeshTest.verifyMesh(testCase,graded,mesh,3);
        end

        function testGradedTwoElements(testCase)
            % the grading point [0 0] is no element mark on a mesh with
            % two elements
            mesh = MeshTest.domainMesh(struct('name','Square'));
            testCase.assertEqual(mesh.nrElems,2);
            H = 0.25;

            graded = mesh.refine('graded',0.5,H,[0 0],1/2);
            beta = mesh.refine('gradedbeta',2/3,[0 0]);

            length4ed = graded.length4ed;
            testCase.verifyLessThanOrEqual(max(length4ed(graded.ed4e),[],2),H + testCase.tol);
            testCase.verifyNumElements(graded.parent4e,graded.nrElems);
            testCase.verifyEqual(beta.nrElems,8);
            MeshTest.verifyMesh(testCase,graded,mesh,1);
            MeshTest.verifyMesh(testCase,beta,mesh,1);
        end

        function testGradedPoints(testCase)
            % several grading points with one exponent and radius each
            % or common ones
            mesh = MeshTest.domainMesh(struct('name','Lshape'));
            points = [0 0;-1 1];
            H = 0.25;

            common = mesh.refine('graded',0.5,H,points);
            each = mesh.refine('graded',[0.5;0.7],H,points,[1;1/2]);
            red = mesh.refine('gradedred',[0.25 0.4],points);

            length4ed = common.length4ed;
            hT = max(length4ed(common.ed4e),[],2);
            for k = 1 : size(points,1)
                rT = sqrt(sum((common.point4e([1 1]/3) - points(k,:)).^2,2));
                testCase.verifyLessThanOrEqual(hT,H*rT.^(1-0.5) + testCase.tol,sprintf('point %d',k));
            end
            testCase.verifyLessThan(each.nrElems,common.nrElems,'weaker grading');
            % the new nodes at the second point divide its edges by 0.4
            corner = find(all(mesh.c4n == points(2,:),2));
            atCorner = find(any(mesh.n4ed == corner,2));
            other = sum(mesh.n4ed(atCorner,:),2) - corner;
            expected = points(2,:) + 0.4*(mesh.c4n(other,:) - points(2,:));
            testCase.verifyEqual(red.c4n(mesh.nrNodes + atCorner,:),expected,'AbsTol',testCase.tol);
            MeshTest.verifyMesh(testCase,common,mesh,3);
            MeshTest.verifyMesh(testCase,each,mesh,3);
            MeshTest.verifyMesh(testCase,red,mesh,3);
        end

        function testSaveLoadExact(testCase)
            % curved boundaries give coordinates with all digits
            mesh = Mesh(Geometry(struct('name','CircularInclusion'))).refine('uniform');
            file = [tempname,'_mesh'];
            cleanup = onCleanup(@()delete([file,'_*.dat']));

            mesh.save(file);
            loaded = Mesh(file);

            testCase.verifyEqual(loaded.c4n,mesh.c4n,'coordinates');
            testCase.verifyEqual(loaded.n4e,mesh.n4e,'elements');
            testCase.verifyEqual(loaded.n4Sb,mesh.n4Sb,'symmetry pairs');
        end

        function testSymmetryAxisRequired(testCase)
            % files with node pairs only, without the fixed axis
            mesh = Mesh(Geometry(struct('name','Square','bdry','SDSN'))).refine('uniform');
            file = [tempname,'_mesh'];
            cleanup = onCleanup(@()delete([file,'_*.dat']));
            mesh.save(file);
            writematrix(mesh.n4Sb(:,1:2),[file,'_n4Sb.dat'],'Delimiter',' ');

            testCase.verifyError(@()Mesh(file),'fem4m:Mesh:invalidN4Sb','Mesh(NAME)');
            testCase.verifyError(@()Geometry(struct('name',file)),'fem4m:Geometry:invalidN4Sb','Geometry');
            testCase.verifyError(@()MeshPoly(file),'fem4m:MeshPoly:invalidN4Sb','MeshPoly(NAME)');
        end

        function testFarSymmetryAxis3d(testCase)
            % far-off coordinates turned by a quarter keep their axes
            mesh = Mesh3d(Geometry(struct('name','Cube','bdry','SSSSSS', ...
                'midpoint',[0 1000 0.5]))).refine('uniform');
            file = [tempname,'_mesh'];
            cleanup = onCleanup(@()delete([file,'_*.dat']));
            mesh.save(file);

            turned = Geometry(struct('name',file,'rotation',pi/2));
            triples = Mesh3d(file);

            % x faces turn into y faces and y faces into x faces
            swapped = [0 2 1 3]';
            expected = swapped(double(mesh.n4Sb(:,4))+1);
            testCase.verifyEqual(turned.n4Sb(:,4),expected,'rotated file');
            testCase.verifyEqual(triples.n4Sb,mesh.n4Sb,'reloaded');
            % face triples without the fixed axis
            writematrix(mesh.n4Sb(:,1:3),[file,'_n4Sb.dat'],'Delimiter',' ');
            testCase.verifyError(@()Mesh3d(file),'fem4m:Mesh3d:invalidN4Sb','Mesh3d(NAME)');
            testCase.verifyError(@()Geometry(struct('name',file,'rotation',pi/2)), ...
                'fem4m:Geometry:invalidN4Sb','rotated Geometry');
        end

        function testColoring(testCase,domain)
            mesh = MeshTest.domainMesh(domain).refine('uniform');

            color4n = mesh.color4n;
            color4n2 = mesh.color4n(2);
            color4e = mesh.color4e;

            testCase.verifyGreaterThanOrEqual(color4n,1);
            testCase.verifyTrue(all(color4n(mesh.n4ed(:,1)) ~= color4n(mesh.n4ed(:,2))), ...
                'nodes of an edge');
            % with two levels, nodes two edges apart also differ
            A = spones(mesh.ed4n);
            [i,j] = find(A + A*A);
            far = i ~= j;
            testCase.verifyTrue(all(color4n2(i(far)) ~= color4n2(j(far))),'two levels');
            patch4n = mesh.patch4n;
            for k = 1 : mesh.nrNodes
                testCase.verifyNumElements(unique(color4e(patch4n{k})),numel(patch4n{k}), ...
                    sprintf('elements at node %d',k));
            end
        end

        function testPatches(testCase)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));

            patch4e = mesh.patch4e;
            patch4n = mesh.patch4n;

            testCase.verifySize(patch4e,[mesh.nrElems,1]);
            testCase.verifySize(patch4n,[mesh.nrNodes,1]);
            for k = 1 : mesh.nrNodes
                testCase.verifyEqual(patch4n{k}(:),find(any(mesh.n4e == k,2)), ...
                    sprintf('patch4n{%d}',k));
            end
            for T = 1 : mesh.nrElems
                expected = find(any(ismember(mesh.n4e,mesh.n4e(T,:)),2));
                testCase.verifyEqual(patch4e{T}(:),expected,sprintf('patch4e{%d}',T));
            end
        end

        function testRemoveNode(testCase)
            mesh = MeshTest.domainMesh(struct('name','Square','bdry','DNDN')).refine('uniform',2);
            area4e = mesh.area4e;
            % an interior node and a node inside a Neumann side
            interior = find(all(abs(mesh.c4n - 0.5) < testCase.tol,2));
            neumann = setdiff(unique(mesh.n4Nb),unique(mesh.n4Db));
            nodes = [interior,neumann(1)];
            patch = find(any(ismember(mesh.n4e,nodes),2));

            removed = mesh.removeNode(nodes);

            testCase.verifyEqual(removed.nrNodes,mesh.nrNodes - 2);
            testCase.verifyEqual(removed.nrElems,mesh.nrElems - numel(patch));
            testCase.verifyGreaterThan(size(removed.n4Nb,1),size(mesh.n4Nb,1), ...
                'exposed sides at a Neumann node are Neumann');
            MeshTest.verifyMesh(testCase,removed,removed,1 - sum(area4e(patch)));
        end

        function testRemoveEdgeAddTriangle(testCase)
            mesh = MeshTest.domainMesh(struct('name','Square')).refine('uniform');
            area4e = mesh.area4e;
            edge = mesh.ed4Db(1);
            elem = mesh.e4ed(edge,1);

            removed = mesh.removeEdge(edge);
            restored = removed.addTriangle(mesh.n4e(elem,:));

            testCase.verifyEqual(removed.nrElems,mesh.nrElems - 1);
            MeshTest.verifyMesh(testCase,removed,removed,1 - area4e(elem));
            testCase.verifyEqual(restored.nrElems,mesh.nrElems);
            testCase.verifyEqual(restored.nrEdges,mesh.nrEdges);
            MeshTest.verifyMesh(testCase,restored,mesh,1);
        end

        function testRemoveSymmetryEdge(testCase)
            mesh = MeshTest.domainMesh(struct('name','Square','bdry','SSSS')).refine('uniform');
            area4e = mesh.area4e;
            edge = mesh.ed4Sb(1);
            elem = mesh.e4ed(edge,1);

            removed = mesh.removeEdge(edge);

            testCase.verifyEqual(removed.nrElems,mesh.nrElems - 1);
            testCase.verifyEqual(size(removed.n4Sb),size(mesh.n4Sb) + [1 0], ...
                'exposed sides at a symmetry edge are symmetry edges');
            testCase.verifyEqual(removed.n4Sb(end-1:end,3),mesh.n4Sb([1 1],3), ...
                'exposed sides inherit the symmetry axis');
            MeshTest.verifyMesh(testCase,removed,removed,1 - area4e(elem));
        end

        function testConnectivity3d(testCase)
            mesh = UnitTestCase.cubeMesh();

            e4n = mesh.e4n;
            f4n = mesh.f4n;
            ed4n = mesh.ed4n;
            e4ed = mesh.e4ed;
            f4ed = mesh.f4ed;

            testCase.verifySize(e4n,[mesh.nrNodes,mesh.nrElems]);
            testCase.verifyEqual(nnz(e4n),4*mesh.nrElems);
            testCase.verifyEqual(nnz(f4n),3*mesh.nrFaces);
            testCase.verifyEqual(nnz(ed4n),2*mesh.nrEdges);
            testCase.verifyEqual(nnz(e4ed),6*mesh.nrElems);
            testCase.verifyEqual(nnz(f4ed),3*mesh.nrFaces);
            testCase.verifyEqual(full(ed4n - ed4n.'),zeros(mesh.nrNodes),'ed4n symmetric');
            for j = 1 : 4
                testCase.verifyEqual(full(e4n(sub2ind(size(e4n),double(mesh.n4e(:,j)), ...
                    (1:mesh.nrElems)'))),(1:mesh.nrElems)',sprintf('e4n of n4e(:,%d)',j));
            end
            testCase.verifyEqual(full(ed4n(sub2ind(size(ed4n),double(mesh.n4ed(:,1)), ...
                double(mesh.n4ed(:,2))))),(1:mesh.nrEdges)','ed4n of n4ed');
            for j = 1 : 6
                testCase.verifyEqual(full(e4ed(sub2ind(size(e4ed),double(mesh.ed4e(:,j)), ...
                    (1:mesh.nrElems)'))),(1:mesh.nrElems)',sprintf('e4ed of ed4e(:,%d)',j));
            end
            for j = 1 : 3
                testCase.verifyEqual(full(f4ed(sub2ind(size(f4ed),double(mesh.ed4f(:,j)), ...
                    (1:mesh.nrFaces)'))),(1:mesh.nrFaces)',sprintf('f4ed of ed4f(:,%d)',j));
                testCase.verifyEqual(full(f4n(sub2ind(size(f4n),double(mesh.n4f(:,j)), ...
                    (1:mesh.nrFaces)'))),(1:mesh.nrFaces)',sprintf('f4n of n4f(:,%d)',j));
            end
        end

        function testConnectivity2d(testCase)
            % row n of e4n lists the elements at the node n, as in 3d
            meshes = {Mesh(Geometry(struct('name','Lshape'))).refine('uniform'), ...
                MeshPolyH(Geometry(struct('name','Square','type','H')))};
            for k = 1 : numel(meshes)
                mesh = meshes{k};
                if iscell(mesh.n4e)
                    nodes = [mesh.n4e{:}]';
                    elems = repelem((1:mesh.nrElems)',mesh.nrN4e);
                else
                    nodes = mesh.n4e(:);
                    elems = repmat((1:mesh.nrElems)',size(mesh.n4e,2),1);
                end
                testCase.verifySize(mesh.e4n,[mesh.nrNodes,mesh.nrElems],class(mesh));
                testCase.verifyEqual(nnz(mesh.e4n),numel(nodes),class(mesh));
                testCase.verifyEqual(full(mesh.e4n(sub2ind(size(mesh.e4n),nodes,elems))), ...
                    elems,class(mesh));
            end
        end

        function testConnectivitySizes(testCase)
            % an unused last node keeps the sizes of the node matrices
            mesh = Mesh([0 0;1 0;0 1;5 5],[1 2 3],[1 2;2 3;3 1],[]);
            testCase.verifySize(mesh.e4n,[4 1],'2d e4n');
            testCase.verifySize(mesh.ed4n,[4 4],'2d ed4n');
            mesh = Mesh3d([0 0 0;1 0 0;0 1 0;0 0 1;9 9 9],[1 2 3 4], ...
                [1 3 2;1 2 4;1 4 3;2 3 4],[],[],'t4e',0);
            testCase.verifySize(mesh.e4n,[5 1],'3d e4n');
            testCase.verifySize(mesh.ed4n,[5 5],'3d ed4n');
            testCase.verifySize(mesh.f4n,[5 4],'3d f4n');
            testCase.verifySize(mesh.e4ed,[6 1],'3d e4ed');
            testCase.verifySize(mesh.f4ed,[6 4],'3d f4ed');
            % the 3d incidences are computed on access and store nothing
            props = ?Mesh3d;
            props = props.PropertyList;
            for name = {'e4n','ed4n','e4ed','f4n','f4ed'}
                testCase.verifyTrue(props(strcmp({props.Name},name{1})).Dependent,name{1});
            end
        end

        function testNormalsTangents3d(testCase)
            mesh = UnitTestCase.cubeMesh();
            tol = testCase.tol;
            perm = [1 2 3;1 4 2;2 4 3;1 3 4];

            normal4e = mesh.normal4e;
            tangent4ed = mesh.tangent4ed;

            testCase.verifySize(normal4e,[mesh.nrElems,4,3]);
            mid4e = mesh.point4e([1 1 1]/4);
            weighted = zeros(mesh.nrElems,3);
            for j = 1 : 4
                n = reshape(normal4e(:,j,:),[],3);
                P = cell(1,3);
                for i = 1 : 3
                    P{i} = mesh.c4n(mesh.n4e(:,perm(j,i)),:);
                end
                testCase.verifyEqual(sum(n.^2,2),ones(mesh.nrElems,1),'AbsTol',tol, ...
                    sprintf('unit normal %d',j));
                testCase.verifyEqual(sum(n.*(P{2}-P{1}),2),zeros(mesh.nrElems,1),'AbsTol',tol, ...
                    sprintf('normal %d orthogonal to the face',j));
                testCase.verifyEqual(sum(n.*(P{3}-P{1}),2),zeros(mesh.nrElems,1),'AbsTol',tol, ...
                    sprintf('normal %d orthogonal to the face',j));
                testCase.verifyGreaterThan(sum(n.*((P{1}+P{2}+P{3})/3 - mid4e),2),0, ...
                    sprintf('normal %d points outward',j));
                areaF = sqrt(sum(cross(P{2}-P{1},P{3}-P{1}).^2,2))/2;
                weighted = weighted + areaF.*n;
            end
            testCase.verifyEqual(weighted,zeros(mesh.nrElems,3),'AbsTol',tol,'sum |F| n_F');
            testCase.verifyEqual(sum(tangent4ed.^2,2),ones(mesh.nrEdges,1),'AbsTol',tol);
            testCase.verifyEqual(tangent4ed.*mesh.length4ed, ...
                mesh.c4n(mesh.n4ed(:,2),:) - mesh.c4n(mesh.n4ed(:,1),:),'AbsTol',tol);
        end

        function testEnumeration3d(testCase)
            meshes = {UnitTestCase.cubeMesh('DNDNDN'), ...
                Mesh3d(Geometry(struct('name','Cube','reftype','crisscross','bdry','DNDNDN'))), ...
                Mesh3d(Geometry(struct('name','Lshape3d'))), ...
                UnitTestCase.referenceTetrahedron()};

            for k = 1 : numel(meshes)
                MeshTest.verifyEnumeration3d(testCase,meshes{k});
            end
        end

        function testUnlistedBoundaryFaces3d(testCase)
            mesh = UnitTestCase.cubeMesh('DNDNDN');

            % without the Neumann triples their faces are in no list
            testCase.verifyError(@()Mesh3d(mesh.c4n,mesh.n4e,mesh.n4Db,[],[],'t4e',mesh.t4e), ...
                'fem4m:Mesh3d:unmatchedFaces');
        end

        function testBisectionLabels3d(testCase)
            % the stored labels of the 3d geometries keep every uniform
            % sweep conforming: one sweep doubles the elements
            names = {'Kuhn cube','crisscross cube','Cube files','Lshape3d','Fichera'};
            meshes = {Mesh3d(Geometry(struct('name','Cube'))), ...
                Mesh3d(Geometry(struct('name','Cube','reftype','crisscross'))), ...
                Mesh3d('Cube'), ...
                Mesh3d(Geometry(struct('name','Lshape3d'))), ...
                Mesh3d(Geometry(struct('name','Fichera')))};
            for k = 1 : numel(meshes)
                mesh = meshes{k};
                testCase.verifyTrue(all(mesh.t4e >= 0 & mesh.t4e <= 5),'labels 0..5');
                refined = mesh;
                for sweep = 1 : 3
                    refined = refined.refine(true(refined.nrElems,1));
                    testCase.verifyEqual(refined.nrElems,2^sweep*mesh.nrElems,names{k});
                end
                MeshTest.verifyMesh3d(testCase,refined,mesh,sum(mesh.area4e));
            end
        end

        function testLongestEdgeLabels3d(testCase)
            % the stored labels of the 3d geometries are the longest edge
            % labels
            geoms = {struct('name','Cube'),struct('name','Cube','reftype','crisscross'), ...
                struct('name','Lshape3d'),struct('name','Fichera')};
            for k = 1 : numel(geoms)
                g = Geometry(geoms{k});
                stored = Mesh3d(g);
                computed = Mesh3d(g.c4n,g.n4e,g.n4Db,g.n4Nb,g.n4Sb);
                testCase.verifyEqual(MeshTest.refineSequence(computed), ...
                    MeshTest.refineSequence(stored),'AbsTol',1e-12,geoms{k}.name);
            end
            % a Delaunay mesh with all four types
            c4n = [0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1; ...
                mod((1:40)'*[0.618034 0.414214 0.732051],1)];
            n4e = delaunay(c4n);
            mesh = Mesh3d(c4n,n4e,Mesh3d.getBoundary(n4e),[]);
            kinds = {'A',[0 3];'P',[1 4];'M',[6 8];'O',[7 9]};
            for k = 1 : size(kinds,1)
                testCase.verifyTrue(any(ismember(mesh.t4e,kinds{k,2})),kinds{k,1});
            end
            length4ed = mesh.length4ed;
            testCase.verifyEqual(vecnorm(mesh.c4n(mesh.n4e(:,4),:) - mesh.c4n(mesh.n4e(:,1),:),2,2), ...
                max(length4ed(mesh.ed4e),[],2),'RelTol',1e-12,'longest edge first');
            uniform = mesh.refine('uniform',2);
            testCase.verifyEqual(accumarray(double(uniform.parent4e(:)),1), ...
                64*ones(mesh.nrElems,1),'two uniform rounds');
            MeshTest.verifyMesh3d(testCase,uniform,mesh,1);
            refined = mesh;
            for k = 1 : 4
                refined = refined.refine(MeshTest.hashMarks(refined));
            end
            MeshTest.verifyMesh3d(testCase,refined,mesh,1);
            % one element of every type alone, and two elements with
            % equally long edges: all elements of their refinements are in
            % the bisection tree of Arnold, Mukherjee, Pouly
            singles = cell(size(kinds,1),3);
            for k = 1 : size(kinds,1)
                e = find(ismember(mesh.t4e,kinds{k,2}),1);
                singles(k,:) = {mesh.c4n(mesh.n4e(e,:),:),kinds{k,1},kinds{k,2}};
            end
            singles(end+1,:) = {[0 0 0;1 0 0;0 1 0;0 0 1],'P with equal edges',[1 4]};
            singles(end+1,:) = {[0 0 0;2 0 0;1 1 0;1 0 1],'O with equal edges',[7 9]};
            faces = [1 2 3;1 2 4;1 3 4;2 3 4];
            for k = 1 : size(singles,1)
                [X,name,labels] = singles{k,:};
                single = Mesh3d(X,[1 2 3 4],faces,[]);
                testCase.verifyTrue(ismember(single.t4e,labels),name);
                tree = MeshTest.ampTree(X,9);
                for sweep = 1 : 3
                    single = single.refine(true(single.nrElems,1));
                    testCase.verifyTrue(all(ismember(MeshTest.tetRows(single.c4n,single.n4e), ...
                        tree,'rows')),sprintf('%s, sweep %d',name,sweep));
                end
            end
        end

        function testInconsistentLabels3d(testCase)
            % labels from the coloring with the wrong type: the closure
            % cannot see the face mismatch, the enumeration of the
            % refined mesh does
            g = Geometry(struct('name','cube','reftype','crisscross'));
            colored = Mesh3d(g.c4n,g.n4e,g.n4Db,g.n4Nb,'init','coloring');
            bad = Mesh3d(colored.c4n,colored.n4e,colored.n4Db,colored.n4Nb,[], ...
                't4e',ones(colored.nrElems,1));

            function refineThrice()
                m = bad;
                for k = 1 : 3
                    m = m.refine(true(m.nrElems,1));
                end
            end

            testCase.verifyError(@refineThrice,'fem4m:Mesh3d:unmatchedFaces');
        end

        function testLabelRoundTrips3d(testCase)
            mesh = Mesh3d(Geometry(struct('name','Lshape3d','bdry','D')));
            mesh = mesh.refine(MeshTest.hashMarks(mesh));
            % the same mesh given in bisection order, as stored, mirrored
            % and saved
            maubach = double(mesh.n4e);
            flip = mesh.t4e >= 3;
            maubach(flip,[2 3]) = maubach(flip,[3 2]);
            file = [tempname,'_mesh'];
            mesh.save(file);
            mirrored = mesh.c4n.*[-1 1 1];
            variants = {
                Mesh3d(mesh.c4n,maubach,mesh.n4Db,mesh.n4Nb,[],'t4e',mod(mesh.t4e,3)), [1 1 1]
                Mesh3d(mesh.c4n,mesh.n4e,mesh.n4Db,mesh.n4Nb,[],'t4e',mesh.t4e), [1 1 1]
                Mesh3d(mirrored,mesh.n4e,mesh.n4Db,mesh.n4Nb,[],'t4e',mesh.t4e), [-1 1 1]
                Mesh3d(file), [1 1 1]};
            delete([file,'_*.dat']);
            reference = MeshTest.refineSequence(mesh);
            for k = 1 : size(variants,1)
                testCase.verifyGreaterThan(variants{k,1}.area4e,0,'positive orientation');
                testCase.verifyEqual(MeshTest.refineSequence(variants{k,1},variants{k,2}), ...
                    reference,'AbsTol',1e-12,sprintf('variant %d',k));
            end
        end

        function testSaveLoad3d(testCase)
            % labels and symmetry axes survive save and load, and empty
            % lists delete the files of an earlier save
            mesh = Mesh3d(Geometry(struct('name','Cube','bdry','DSNSDS'))).refine(true(6,1));
            dirichlet = Mesh3d(Geometry(struct('name','Cube','bdry','DDDDDD')));
            file = [tempname,'_mesh'];
            cleanup = onCleanup(@()delete([file,'_*.dat']));

            mesh.save(file);
            loaded = Mesh3d(file);
            dirichlet.save(file);
            reloaded = Mesh3d(file);

            testCase.verifyEqual(loaded.n4e,mesh.n4e,'elements');
            testCase.verifyEqual(loaded.t4e,mesh.t4e,'labels');
            testCase.verifyEqual(loaded.n4Sb,mesh.n4Sb,'n4Sb with axis');
            testCase.verifyEqual(reloaded.n4Db,dirichlet.n4Db);
            testCase.verifyEmpty(reloaded.n4Nb,'no stale Neumann faces');
            testCase.verifyEmpty(reloaded.n4Sb,'no stale symmetry faces');
        end

        function testColoringLabels3d(testCase)
            % an unstructured mesh with labels from the vertex coloring,
            % also with permuted node ids
            c4n = [0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1; ...
                0.4 0.5 0.6;0.7 0.2 0.3;0.2 0.8 0.4;0.6 0.6 0.2;0.5 0.3 0.8];
            n4e = delaunay(c4n);
            perm = [5 12 1 9 3 13 7 2 11 6 10 4 8];
            inverse(perm) = 1:numel(perm);
            permuted = inverse(n4e);
            meshes = {Mesh3d(c4n,n4e,Mesh3d.getBoundary(n4e),[],'init','coloring'), ...
                Mesh3d(c4n(perm,:),permuted,Mesh3d.getBoundary(permuted),[],'init','coloring')};
            for k = 1 : numel(meshes)
                testCase.verifyTrue(all(meshes{k}.t4e == 0 | meshes{k}.t4e == 3),'type 0');
                refined = meshes{k};
                for level = 1 : 4
                    refined = refined.refine(MeshTest.hashMarks(refined));
                end
                MeshTest.verifyMesh3d(testCase,refined,meshes{k},1);
            end
        end

        function testEdgeMarks3d(testCase)
            mesh = UnitTestCase.cubeMesh('DNDNDN');
            marked4ed = false(mesh.nrEdges,1);
            marked4ed(mesh.ed4e(1,2)) = true;   % not the refinement edge

            refined = mesh.refine(marked4ed);

            midpoint = mean(mesh.c4n(mesh.n4ed(mesh.ed4e(1,2),:),:),1);
            testCase.verifyEqual(min(vecnorm(refined.c4n - midpoint,2,2)),0, ...
                'AbsTol',testCase.tol,'the edge is halved');
            MeshTest.verifyMesh3d(testCase,refined,mesh,1);
        end

        function testSymmetryBoundary3d(testCase)
            geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5], ...
                'length',1,'bdry','DSNSDS'));
            % the faces x = max, x = min and z = min
            testCase.verifyEqual(geom.n4Sb(:,4),[1;1;1;1;3;3],'axis of the geometry');
            mesh = Mesh3d(geom);
            testCase.assertNotEmpty(mesh.n4Sb);

            refined = mesh.refine('uniform');
            refined = refined.refine(MeshTest.hashMarks(refined));

            MeshTest.verifyMesh3d(testCase,refined,mesh,1);
            % the axis is the direction of the face normal and inherited
            c1 = refined.c4n(refined.n4Sb(:,1),:);
            normal = cross(refined.c4n(refined.n4Sb(:,2),:)-c1,refined.c4n(refined.n4Sb(:,3),:)-c1,2);
            [~,direction] = max(abs(normal),[],2);
            testCase.verifyEqual(double(refined.n4Sb(:,4)),direction,'axis = normal direction');
            area4f = refined.area4f;
            initialArea4f = mesh.area4f;
            testCase.verifyEqual(sum(area4f(refined.f4Sb)),sum(initialArea4f(mesh.f4Sb)), ...
                'AbsTol',testCase.tol,'symmetry boundary area');
        end

        function testConstructorPaths3d(testCase)
            mesh = Mesh3d(Geometry(struct('name','Cube','bdry','DSNSDS')));

            fromStruct = Mesh3d(struct('name','Cube','bdry','DSNSDS'));
            explicit = Mesh3d(mesh.c4n,mesh.n4e,mesh.n4Db,mesh.n4Nb,mesh.n4Sb,'t4e',mesh.t4e);

            testCase.verifyEqual(fromStruct.n4e,mesh.n4e,'options struct');
            testCase.verifyEqual(fromStruct.n4Sb,mesh.n4Sb,'options struct');
            testCase.verifyEqual(explicit.n4e,mesh.n4e,'explicit arrays');
            testCase.verifyEqual(explicit.n4Sb,mesh.n4Sb,'explicit arrays');
            testCase.verifyError(@()Mesh3d(mesh.c4n,mesh.n4e,mesh.n4Db,mesh.n4Nb, ...
                mesh.n4Sb(:,1:3),'t4e',mesh.t4e),'fem4m:Mesh3d:invalidN4Sb','triples without the axis');
            text = evalc('disp(explicit)');
            testCase.verifyTrue(contains(text,'6 elements') && contains(text,'6 symmetry'),'disp');
            color4n = mesh.color4n;
            testCase.verifyTrue(all(color4n(mesh.n4ed(:,1)) ~= color4n(mesh.n4ed(:,2))),'coloring');
        end

        function testLargeNodeIds3d(testCase)
            % with more than 2^(64/3) nodes the face keys overflow uint64
            % and the enumeration sorts rows instead
            mesh = UnitTestCase.cubeMesh();
            nrNodes = 3e6;
            ids = nrNodes - mesh.nrNodes + (1:mesh.nrNodes)';
            c4n = zeros(nrNodes,3);
            c4n(ids,:) = mesh.c4n;

            large = Mesh3d(c4n,ids(mesh.n4e),ids(mesh.n4Db),ids(mesh.n4Nb),[],'t4e',mesh.t4e);

            testCase.verifyEqual(double(large.n4f),ids(mesh.n4f),'faces');
            testCase.verifyEqual(double(large.n4ed),ids(mesh.n4ed),'edges');
            testCase.verifyEqual(large.f4e,mesh.f4e,'faces of the elements');
            testCase.verifyEqual(large.ed4f,mesh.ed4f,'edges of the faces');
            testCase.verifyEqual(large.f4Db,mesh.f4Db,'Dirichlet faces');
        end

        function testRotatedSymmetryBoundary3d(testCase)
            opts = struct('name','Cube','bdry','DSNSDS','rotation',pi/2);
            geom = Geometry(opts);
            % the x faces turn into y faces
            testCase.verifyEqual(geom.n4Sb(:,4),[2;2;2;2;3;3],'quarter turn');

            opts.rotation = pi/4;
            geom = Geometry(opts);
            mesh = Mesh3d(geom);

            % general normals get axis 0
            testCase.verifyEqual(geom.n4Sb(:,4),[0;0;0;0;3;3],'eighth turn');
            testCase.verifyEqual(sort(double(mesh.n4Sb(:,4))),[0;0;0;0;3;3],'mesh axis');
            refined = mesh.refine('uniform');
            testCase.verifyEqual(accumarray(double(refined.n4Sb(:,4))+1,1),[16;0;0;8], ...
                'refined axis');
        end

        function testProlongationData3d(testCase)
            mesh = UnitTestCase.cubeMesh();

            refined = mesh.refine(MeshTest.hashMarks(mesh));

            % the P1 prolongation interpolates the coordinates and every
            % element lies in its ancestor
            testCase.verifySize(refined.prolongation,[refined.nrNodes,mesh.nrNodes]);
            testCase.verifyEqual(refined.prolongation*mesh.c4n,refined.c4n,'AbsTol',testCase.tol);
            testCase.verifyEqual(accumarray(double(refined.parent4e(:)),refined.area4e), ...
                mesh.area4e,'AbsTol',testCase.tol);
            uniform = mesh.refine('uniform',2);
            testCase.verifyEqual(uniform.prolongation*mesh.c4n,uniform.c4n,'AbsTol',testCase.tol);
            testCase.verifyEqual(accumarray(double(uniform.parent4e(:)),uniform.area4e), ...
                mesh.area4e,'AbsTol',testCase.tol);
        end

        function testBoundaryAndDeleteNodes3d(testCase)
            mesh = UnitTestCase.cubeMesh();
            node = 1;

            [n4f,e4f] = Mesh3d.getBoundary(mesh.n4e);
            [c4n,n4e,t4e,n4bd] = Mesh3d.deleteNodes(mesh.c4n,mesh.n4e,mesh.t4e,node);

            boundary = mesh.e4f(:,2) == 0;
            testCase.verifyEqual(size(n4f,1),nnz(boundary));
            testCase.verifyEqual(sortrows(double(n4f)),sortrows(sort(double(mesh.n4f(boundary,:)),2)));
            testCase.verifyEqual(sort(double(e4f)),sort(double(mesh.e4f(boundary,1))));
            kept = ~any(mesh.n4e == node,2);
            testCase.verifyEqual(c4n,mesh.c4n([1:node-1,node+1:end],:));
            testCase.verifyEqual(size(n4e,1),nnz(kept));
            testCase.verifyEqual(numel(t4e),nnz(kept));
            testCase.verifyEqual(c4n(n4e,:),mesh.c4n(mesh.n4e(kept,:),:),'node coordinates');
            testCase.verifyEqual(n4bd,Mesh3d.getBoundary(n4e));
        end

    end

    methods (Static)

        function mesh = domainMesh(domain)
            %domainMesh - Initial mesh of the domain without refinement
            opts = rmfield(domain,intersect(fieldnames(domain),{'area'}));
            mesh = Mesh(Geometry(opts));
        end

        function verifyMesh(testCase,mesh,initialMesh,area)
            %verifyMesh - Verify invariants and consistency of a 2d mesh
            %   MeshTest.verifyMesh(TESTCASE,MESH,INITIALMESH,AREA)
            %   verifies that MESH is a conforming, positively oriented
            %   triangulation of the domain of area AREA with consistent
            %   edge and boundary data that has the same boundary length
            %   and Euler characteristic as INITIALMESH.
            tol = testCase.tol;
            area4e = mesh.area4e;
            testCase.verifyGreaterThan(area4e,0,'positive area');
            testCase.verifyEqual(sum(area4e),area,'AbsTol',tol,'total area');
            testCase.verifyEqual(mesh.nrNodes,size(mesh.c4n,1));
            testCase.verifyEqual(mesh.nrElems,size(mesh.n4e,1));
            testCase.verifyEqual(mesh.nrEdges,size(mesh.n4ed,1));
            testCase.verifyEqual(sort(unique(mesh.n4e(:))),(1:mesh.nrNodes)', ...
                'every node belongs to an element');

            % edges of the elements: edge j connects node j and j+1
            for j = 1 : 3
                edges = full(mesh.ed4n(sub2ind(size(mesh.ed4n), ...
                    mesh.n4e(:,j),mesh.n4e(:,mod(j,3)+1))));
                testCase.verifyEqual(edges,mesh.ed4e(:,j),sprintf('ed4e(:,%d)',j));
            end
            testCase.verifyEqual(full(mesh.ed4n(sub2ind(size(mesh.ed4n), ...
                mesh.n4ed(:,1),mesh.n4ed(:,2)))),(1:mesh.nrEdges)','ed4n');

            % adjacent elements of the edges
            interior = mesh.e4ed(:,2) > 0;
            testCase.verifyGreaterThan(mesh.e4ed(:,1),0,'e4ed');
            testCase.verifyEqual(mesh.ed4e(sub2ind(size(mesh.ed4e), ...
                mesh.e4ed(:,1),mesh.edNr4ed(:,1))),(1:mesh.nrEdges)','e4ed(:,1)');
            testCase.verifyEqual(mesh.ed4e(sub2ind(size(mesh.ed4e), ...
                mesh.e4ed(interior,2),mesh.edNr4ed(interior,2))), ...
                find(interior),'e4ed(:,2)');

            % conformity: interior edges belong to two elements, boundary
            % edges to one element and to exactly one boundary part
            testCase.verifyEqual(accumarray(mesh.ed4e(:),1,[mesh.nrEdges,1]), ...
                1 + double(interior),'conforming');
            ed4bd = [mesh.ed4Db;mesh.ed4Nb;mesh.ed4Sb];
            testCase.verifyEqual(sort(ed4bd),find(~interior),'boundary edges');

            % boundary edges are oriented counterclockwise
            n4Sb = mesh.n4Sb;
            if ~isempty(n4Sb)
                n4Sb = n4Sb(:,1:2);
            end
            n4bd = [mesh.n4Db;mesh.n4Nb;n4Sb];
            T = mesh.e4ed(ed4bd,1);
            [~,k] = max(mesh.n4e(T,:) == n4bd(:,1),[],2);
            testCase.verifyEqual(mesh.n4e(sub2ind(size(mesh.n4e),T,mod(k,3)+1)),n4bd(:,2), ...
                'boundary orientation');
            testCase.verifyEqual(mesh.n4ed(ed4bd,:),n4bd(:,1:2),'n4ed of boundary');

            % boundary length and topology of the initial mesh
            length4ed = mesh.length4ed;
            initialLength4ed = initialMesh.length4ed;
            testCase.verifyEqual(sum(length4ed(~interior)), ...
                sum(initialLength4ed(initialMesh.e4ed(:,2)==0)),'AbsTol',tol, ...
                'boundary length');
            testCase.verifyEqual(sum(length4ed(mesh.ed4Db)), ...
                sum(initialLength4ed(initialMesh.ed4Db)),'AbsTol',tol, ...
                'Dirichlet boundary length');
            testCase.verifyEqual(mesh.nrNodes - mesh.nrEdges + mesh.nrElems, ...
                initialMesh.nrNodes - initialMesh.nrEdges + initialMesh.nrElems, ...
                'Euler characteristic');
        end

        function marked4e = hashMarks(mesh)
            %hashMarks - Pseudo-random element marks from the centroids
            %   The marks do not depend on the element order: about a
            %   quarter of the elements, chosen by a hash of the centroid.
            h = abs(sin(mesh.mid4e*[12.9898;78.233;37.719])*43758.5453);
            marked4e = mod(floor(1000*(h - floor(h))),4) == 0;
            if ~any(marked4e)
                marked4e(1) = true;
            end
        end

        function rows = tetRows(c4n,n4e)
            %tetRows - Sorted vertex coordinates of tetrahedra
            %   ROWS(e,:) lists the vertices of element e sorted by their
            %   coordinates, rounded to 10 digits.
            X = c4n(double(n4e)',:);
            rows = zeros(size(n4e,1),12);
            for e = 1 : size(n4e,1)
                rows(e,:) = reshape(sortrows(X(4*e-3:4*e,:))',1,[]);
            end
            rows = round(rows*1e10)/1e10;
        end

        function rows = ampTree(X,nGen)
            %ampTree - Bisection tree of Arnold, Mukherjee, Pouly
            %   ROWS = ampTree(X,NGEN) lists (as tetRows) the tetrahedra
            %   of the generations 0..NGEN of the bisection tree of the
            %   tetrahedron X (4x3), with the longest edges as refinement
            %   edge and marked edges of the faces, equal lengths ordered
            %   by the vertex numbers as in Mesh3d. A tetrahedron has the
            %   vertices X, the refinement edge ref, the marked edge
            %   mark(k,:) of the face opposite vertex k and a flag.
            pairs = nchoosek(1:4,2);
            len = vecnorm(X(pairs(:,1),:) - X(pairs(:,2),:),2,2);
            % the last of equal lengths: the largest sorted vertex pair
            r = find(len == max(len),1,'last');
            tet = struct('X',X,'ref',pairs(r,:),'mark',zeros(4,2),'flag',false);
            for k = 1 : 4
                inFace = find(~any(pairs == k,2));
                j = find(len(inFace) == max(len(inFace)),1,'last');
                tet.mark(k,:) = pairs(inFace(j),:);
            end
            tets = tet;
            rows = MeshTest.tetRows(X,1:4);
            pos = zeros(1,5);
            for g = 1 : nGen
                children = repmat(tet,1,0);
                for t = tets
                    a = t.ref(1);
                    b = t.ref(2);
                    cd = setdiff(1:4,[a b]);
                    ma = t.mark(b,:);
                    mb = t.mark(a,:);
                    % type P: the marks of the faces acd and bcd meet
                    planar = ~all(ismember(ma,cd)) && ~all(ismember(mb,cd)) ...
                        && setdiff(ma,a) == setdiff(mb,b);
                    if planar && t.flag
                        newMark = [5,setdiff(ma,a)];
                    else
                        newMark = cd;
                    end
                    Y = [t.X;(t.X(a,:) + t.X(b,:))/2];
                    for u = [a b]
                        v = [u 5 cd];
                        mark = [newMark;t.mark(a+b-u,:);u cd(2);u cd(1)];
                        pos(v) = 1 : 4;
                        child = struct('X',Y(v,:),'ref',pos(mark(2,:)), ...
                            'mark',pos(mark),'flag',planar && ~t.flag);
                        children(end+1) = child; %#ok<AGROW>
                    end
                end
                tets = children;
                for t = tets
                    rows = [rows;MeshTest.tetRows(t.X,1:4)]; %#ok<AGROW>
                end
            end
        end

        function sets = refineSequence(mesh,scale)
            %refineSequence - Element sets of a fixed refinement sequence
            %   SETS lists, after two marked and one uniform refinement,
            %   the sorted vertex coordinates (scaled by SCALE) of every
            %   element in sorted order, so equal meshes give equal sets.
            if nargin < 2
                scale = [1 1 1];
            end
            for k = 1 : 2
                h = abs(sin((mesh.mid4e.*scale)*[12.9898;78.233;37.719])*43758.5453);
                mesh = mesh.refine(mod(floor(1000*(h - floor(h))),4) == 0);
            end
            mesh = mesh.refine('uniform');
            X = mesh.c4n(mesh.n4e',:).*scale;
            sets = zeros(mesh.nrElems,12);
            for e = 1 : mesh.nrElems
                sets(e,:) = reshape(sortrows(X(4*e-3:4*e,:))',1,[]);
            end
            sets = sortrows(round(sets*1e10)/1e10);
        end

        function verifyEnumeration3d(testCase,mesh)
            %verifyEnumeration3d - Verify the face and edge numbering of a 3d mesh
            %   MeshTest.verifyEnumeration3d(TESTCASE,MESH) verifies that
            %   the faces and edges of MESH are unique with sorted node
            %   ids, that f4e, ed4e, e4f, fNr4f,
            %   ed4f and the boundary ids are consistent with them, and
            %   that the boundary triples point outward.
            localFaces = [1 2 3;1 2 4;4 2 3;1 3 4];
            localEdges = [1 2;2 3;3 1;1 4;2 4;3 4];
            n4e = double(mesh.n4e);
            n4f = double(mesh.n4f);
            n4ed = double(mesh.n4ed);
            testCase.verifyEqual(n4f,sort(n4f,2),'n4f rows sorted');
            testCase.verifyEqual(size(unique(n4f,'rows'),1),mesh.nrFaces,'n4f unique');
            testCase.verifyEqual(n4ed,sort(n4ed,2),'n4ed rows sorted');
            testCase.verifyEqual(size(unique(n4ed,'rows'),1),mesh.nrEdges,'n4ed unique');
            for k = 1 : 4
                testCase.verifyEqual(n4f(mesh.f4e(:,k),:),sort(n4e(:,localFaces(k,:)),2), ...
                    sprintf('f4e(:,%d)',k));
            end
            for j = 1 : 6
                testCase.verifyEqual(n4ed(mesh.ed4e(:,j),:),sort(n4e(:,localEdges(j,:)),2), ...
                    sprintf('ed4e(:,%d)',j));
            end
            for j = 1 : 3
                testCase.verifyEqual(n4ed(mesh.ed4f(:,j),:),sort(n4f(:,[j,mod(j,3)+1]),2), ...
                    sprintf('ed4f(:,%d)',j));
            end
            faces = (1:mesh.nrFaces)';
            first = reshape(double(mesh.f4e(sub2ind(size(mesh.f4e),double(mesh.e4f(:,1)), ...
                double(mesh.fNr4f(:,1))))),[],1);
            testCase.verifyEqual(first,faces,'e4f(:,1) and fNr4f(:,1)');
            interior = mesh.e4f(:,2) > 0;
            second = reshape(double(mesh.f4e(sub2ind(size(mesh.f4e),double(mesh.e4f(interior,2)), ...
                double(mesh.fNr4f(interior,2))))),[],1);
            testCase.verifyEqual(second,faces(interior),'e4f(:,2) and fNr4f(:,2)');
            testCase.verifyTrue(all(mesh.e4f(interior,1) < mesh.e4f(interior,2)), ...
                'e4f(:,1) is the smaller element');
            lists = {mesh.n4Db,mesh.n4Nb,mesh.n4Sb;mesh.f4Db,mesh.f4Nb,mesh.f4Sb; ...
                mesh.ed4Db,mesh.ed4Nb,mesh.ed4Sb};
            for b = 1 : size(lists,2)
                n4b = double(lists{1,b});
                if isempty(n4b)
                    continue
                end
                testCase.verifyEqual(n4f(lists{2,b},:),sort(n4b(:,1:3),2),'boundary face ids');
                for j = 1 : 3
                    testCase.verifyEqual(n4ed(lists{3,b}(:,j),:), ...
                        sort(n4b(:,[j,mod(j,3)+1]),2),'boundary edge ids');
                end
            end
            MeshTest.verifyOutward3d(testCase,mesh);
        end

        function verifyOutward3d(testCase,mesh)
            %verifyOutward3d - Verify that the boundary triples point outward
            lists = {mesh.n4Db,mesh.n4Nb,mesh.n4Sb;mesh.f4Db,mesh.f4Nb,mesh.f4Sb};
            for b = 1 : size(lists,2)
                n4b = double(lists{1,b});
                if isempty(n4b)
                    continue
                end
                elems = mesh.e4f(lists{2,b},1);
                opposite = sum(double(mesh.n4e(elems,:)),2) - sum(n4b(:,1:3),2);
                c1 = mesh.c4n(n4b(:,1),:);
                normal = cross(mesh.c4n(n4b(:,2),:)-c1,mesh.c4n(n4b(:,3),:)-c1,2);
                testCase.verifyLessThan(dot(normal,mesh.c4n(opposite,:)-c1,2),0, ...
                    'boundary triples point outward');
            end
        end

        function verifyMesh3d(testCase,mesh,initialMesh,volume)
            %verifyMesh3d - Verify invariants and consistency of a 3d mesh
            %   MeshTest.verifyMesh3d(TESTCASE,MESH,INITIALMESH,VOLUME)
            %   verifies that MESH is a conforming, positively oriented
            %   tetrahedral mesh of the domain of volume VOLUME with
            %   consistent face and boundary data that has the same
            %   boundary area and Euler characteristic as INITIALMESH.
            tol = testCase.tol;
            area4e = mesh.area4e;
            testCase.verifyGreaterThan(area4e,0,'positive volume');
            testCase.verifyEqual(sum(area4e),volume,'AbsTol',tol,'total volume');
            testCase.verifyEqual(mesh.nrElems,size(mesh.n4e,1));
            testCase.verifyEqual(mesh.nrFaces,size(mesh.n4f,1));

            % conformity: interior faces belong to two elements, boundary
            % faces to one element and to exactly one boundary part
            interior = mesh.e4f(:,2) > 0;
            testCase.verifyEqual(accumarray(mesh.f4e(:),1,[mesh.nrFaces,1]), ...
                1 + double(interior),'conforming');
            testCase.verifyEqual(sort(double([mesh.f4Db(:);mesh.f4Nb(:);mesh.f4Sb(:)])), ...
                find(~interior),'boundary faces');
            MeshTest.verifyOutward3d(testCase,mesh);
            for j = 1 : 4
                faces = mesh.f4e(:,j);
                testCase.verifyTrue(all(mesh.e4f(faces,1) == (1:mesh.nrElems)' ...
                    | mesh.e4f(faces,2) == (1:mesh.nrElems)'),sprintf('e4f of f4e(:,%d)',j));
            end

            % boundary area and topology of the initial mesh
            area4f = mesh.area4f;
            initialArea4f = initialMesh.area4f;
            testCase.verifyEqual(sum(area4f(~interior)), ...
                sum(initialArea4f(initialMesh.e4f(:,2)==0)),'AbsTol',tol, ...
                'boundary area');
            testCase.verifyEqual(sum(area4f(mesh.f4Db)), ...
                sum(initialArea4f(initialMesh.f4Db)),'AbsTol',tol, ...
                'Dirichlet boundary area');
            testCase.verifyEqual( ...
                mesh.nrNodes - mesh.nrEdges + mesh.nrFaces - mesh.nrElems, ...
                initialMesh.nrNodes - initialMesh.nrEdges + initialMesh.nrFaces ...
                - initialMesh.nrElems,'Euler characteristic');
        end

    end

end
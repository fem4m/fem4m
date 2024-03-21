classdef MarkTest < UnitTestCase
    %MarkTest - Unit tests of the marking strategies for adaptive refinement
    %   The bulk criterion of MarkBulk marks the smallest set of entities
    %   whose squared indicators sum up to theta times the total, MarkMaximum
    %   marks the indicators larger than theta times the maximum, MarkMean
    %   the squared indicators of at least three quarters of their mean per
    %   element, and MarkFixedFraction a fixed number of the largest
    %   indicators. Element indicators also mark the edges of the marked
    %   elements, node indicators the edges and elements at the marked
    %   nodes, in 2d and 3d. MarkUniform and MarkP mark everything for h-
    %   and p-refinement. MarkHP and MarkHP3Edges split the marked
    %   elements into h- and p-refinement by a predictor of the squared
    %   indicators and update the predictor. In an adaptive loop, the
    %   marked edges are refined and the predictor follows the elements to
    %   their children.
    %
    %   See also Mark, MarkBulk, MarkMaximum, MarkMean, MarkFixedFraction,
    %   MarkUniform, MarkP, MarkHP, MarkHP3Edges, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        hpStrategy = struct('MarkHP',@MarkHP,'MarkHP3Edges',@MarkHP3Edges)
        nodeStrategy = struct('MarkBulk',@MarkBulk,'MarkMaximum',@MarkMaximum, ...
            'MarkMean',@MarkMean,'MarkFixedFraction',@MarkFixedFraction)
        theta = struct('small',0.2,'half',0.5,'large',0.9)
    end

    methods (Test)

        function testBulkElementIndicators(testCase,theta)
            fe = MarkTest.feStruct();
            eta4e = MarkTest.indicators(fe.mesh.nrElems);

            mark = MarkBulk(theta).mark(fe,eta4e);

            [~,ind] = sort(eta4e,'descend');
            J = find(cumsum(eta4e(ind).^2) >= theta*sum(eta4e.^2),1,'first');
            expected = zeros(fe.mesh.nrElems,1);
            expected(ind(1:J)) = 1;
            testCase.verifyEqual(mark.marked4e,expected);
            testCase.verifyEqual(mark.marked4ed,MarkTest.edgesOfElements(fe.mesh,expected));
            testCase.verifyEqual(mark.markedP4e,zeros(fe.mesh.nrElems,1));
            % the marked set is minimal: without its smallest member, the
            % bulk criterion fails
            marked = eta4e(expected == 1);
            testCase.verifyGreaterThanOrEqual(sum(marked.^2),theta*sum(eta4e.^2));
            testCase.verifyLessThan(sum(marked.^2) - min(marked)^2,theta*sum(eta4e.^2));
        end

        function testBulkExtremeTheta(testCase)
            fe = MarkTest.feStruct();
            eta4e = MarkTest.indicators(fe.mesh.nrElems);

            none = MarkBulk(0).mark(fe,eta4e);
            everything = MarkBulk(1).mark(fe,eta4e);

            [~,imax] = max(eta4e);
            testCase.verifyEqual(find(none.marked4e),imax,'theta = 0');
            testCase.verifyEqual(everything.marked4e,ones(fe.mesh.nrElems,1),'theta = 1');
            testCase.verifyEqual(everything.marked4ed,ones(fe.mesh.nrEdges,1),'theta = 1');
        end

        function testBulkEdgeIndicators(testCase)
            fe = MarkTest.feStruct();
            eta4ed = MarkTest.indicators(fe.mesh.nrEdges);

            mark = MarkBulk(0.5).mark(fe,eta4ed);

            [~,ind] = sort(eta4ed,'descend');
            J = find(cumsum(eta4ed(ind).^2) >= 0.5*sum(eta4ed.^2),1,'first');
            expected = zeros(fe.mesh.nrEdges,1);
            expected(ind(1:J)) = 1;
            testCase.verifyEqual(mark.marked4ed,expected);
            testCase.verifyEqual(mark.marked4e,zeros(fe.mesh.nrElems,1));
        end

        function testNodeIndicators(testCase,nodeStrategy)
            % a marked node marks all edges and elements at the node, in 2d
            % and in 3d, where refine reads the element marks; with degree
            % nrElems MarkFixedFraction marks a single indicator
            meshes = {PatchTestCase.patchMesh(struct('name','Lshape')), ...
                UnitTestCase.cubeMesh()};
            for k = 1 : numel(meshes)
                mesh = meshes{k};
                fe = struct('mesh',mesh,'degree',mesh.nrElems, ...
                    'degree4e',ones(mesh.nrElems,1));
                eta4n = zeros(mesh.nrNodes,1);
                node = mesh.n4e(1,1);
                eta4n(node) = 1;

                mark = nodeStrategy(0.5).mark(fe,eta4n);

                expectedEdges = zeros(mesh.nrEdges,1);
                expectedEdges(any(mesh.n4ed == node,2)) = 1;
                expectedElems = zeros(mesh.nrElems,1);
                expectedElems(any(mesh.n4e == node,2)) = 1;
                testCase.verifyEqual(mark.marked4ed,expectedEdges,class(mesh));
                testCase.verifyEqual(mark.marked4e,expectedElems,class(mesh));
            end
        end

        function testBulkSeveralIndicators(testCase)
            fe = MarkTest.feStruct();
            eta4e = zeros(fe.mesh.nrElems,1);
            eta4ed = zeros(fe.mesh.nrEdges,1);
            eta4e(2) = 1;
            eta4ed(5) = 1;

            mark = MarkBulk(0.9).mark(fe,eta4e,eta4ed);

            % both largest indicators are needed for 90 percent
            expectedEdges = MarkTest.edgesOfElements(fe.mesh,double((1:fe.mesh.nrElems)' == 2));
            expectedEdges(5) = 1;
            testCase.verifyEqual(find(mark.marked4e),2);
            testCase.verifyEqual(mark.marked4ed,expectedEdges);
        end

        function testMaximum(testCase,theta)
            fe = MarkTest.feStruct();
            eta4e = MarkTest.indicators(fe.mesh.nrElems);

            mark = MarkMaximum(theta).mark(fe,eta4e);

            expected = double(eta4e > theta*max(eta4e));
            testCase.verifyEqual(mark.marked4e,expected);
            testCase.verifyEqual(mark.marked4ed,MarkTest.edgesOfElements(fe.mesh,expected));
        end

        function testMaximumStrictInequality(testCase)
            fe = MarkTest.feStruct();
            eta4e = ones(fe.mesh.nrElems,1);
            eta4e(3) = 2;

            mark = MarkMaximum(0.5).mark(fe,eta4e);

            % indicators equal to theta times the maximum are not marked
            testCase.verifyEqual(find(mark.marked4e),3);
        end

        function testMean(testCase)
            fe = MarkTest.feStruct();
            eta4e = MarkTest.indicators(fe.mesh.nrElems);

            mark = MarkMean().mark(fe,eta4e);

            expected = double(eta4e.^2 >= 0.75*mean(eta4e.^2));
            testCase.verifyEqual(mark.marked4e,expected);
            testCase.verifyEqual(mark.marked4ed,MarkTest.edgesOfElements(fe.mesh,expected));
        end

        function testFixedFraction(testCase)
            for degree = 1 : 2
                fe = MarkTest.feStruct(degree);
                eta4e = MarkTest.indicators(fe.mesh.nrElems);

                mark = MarkFixedFraction().mark(fe,eta4e);

                n = ceil(0.25/degree*fe.mesh.nrElems);
                testCase.verifyEqual(nnz(mark.marked4e),n,sprintf('degree %d',degree));
                testCase.verifyGreaterThanOrEqual(min(eta4e(mark.marked4e == 1)), ...
                    max(eta4e(mark.marked4e == 0)),sprintf('degree %d',degree));
                testCase.verifyEqual(mark.marked4ed, ...
                    MarkTest.edgesOfElements(fe.mesh,mark.marked4e),sprintf('degree %d',degree));
            end
        end

        function testUniformAndP(testCase)
            fe = MarkTest.feStruct();
            eta4e = MarkTest.indicators(fe.mesh.nrElems);
            n = fe.mesh.nrElems;

            uniform = MarkUniform().mark(fe,eta4e);
            p = MarkP().mark(fe,eta4e);

            testCase.verifyEqual(uniform.marked4e,ones(n,1));
            testCase.verifyEqual(uniform.marked4ed,ones(fe.mesh.nrEdges,1));
            testCase.verifyEqual(uniform.markedP4e,zeros(n,1));
            testCase.verifyEqual(p.marked4e,zeros(n,1));
            testCase.verifyEqual(p.marked4ed,zeros(fe.mesh.nrEdges,1));
            testCase.verifyEqual(p.markedP4e,ones(n,1));
        end

        function testHPFirstStepIsP(testCase,hpStrategy)
            fe = MarkTest.feStruct();
            eta4e = MarkTest.indicators(fe.mesh.nrElems);

            mark = hpStrategy();
            mark = mark.mark(fe,eta4e);

            % with the infinite initial predictor, every element in the
            % mean criterion is p-refined
            I = eta4e.^2 >= 0.75*mean(eta4e.^2);
            testCase.verifyEqual(mark.markedP4e,double(I));
            testCase.verifyEqual(mark.marked4e,zeros(fe.mesh.nrElems,1));
            testCase.verifyEqual(mark.marked4ed,zeros(fe.mesh.nrEdges,1));
            testCase.verifyEqual(mark.predict4e(I),mark.gammaP*eta4e(I).^2);
            testCase.verifyEqual(mark.predict4e(~I),inf(nnz(~I),1));
        end

        function testHPSplitAndPredictor(testCase,hpStrategy)
            degree = 2;
            fe = MarkTest.feStruct(degree);
            n = fe.mesh.nrElems;
            eta4e = MarkTest.indicators(n);
            mark = hpStrategy();
            mark.gammaN = 0.9;
            predict4e = 0.5*eta4e.^2;
            predict4e(1:2:end) = 2*eta4e(1:2:end).^2;
            mark.predict4e = predict4e;

            mark = mark.mark(fe,eta4e);

            I = eta4e.^2 >= 0.75*mean(eta4e.^2);
            J = I & eta4e.^2 > predict4e;
            K = I & ~J;
            testCase.assertTrue(any(J) && any(K) && any(~I),'all three cases occur');
            testCase.verifyEqual(mark.marked4e,double(J));
            testCase.verifyEqual(mark.markedP4e,double(K));
            expectedEdges = zeros(fe.mesh.nrEdges,1);
            if isequal(func2str(hpStrategy),'MarkHP')
                expectedEdges(fe.mesh.ed4e(J,1)) = 1;
                hPredict = 1/2*mark.gammaH*(1/2)^degree*eta4e(J).^2;
            else
                expectedEdges(fe.mesh.ed4e(J,:)) = 1;
                hPredict = 1/4*mark.gammaH*(1/2)^(2*degree)*eta4e(J).^2;
            end
            testCase.verifyEqual(mark.marked4ed,expectedEdges);
            testCase.verifyEqual(mark.predict4e(~I),0.9*predict4e(~I),'AbsTol',testCase.tol);
            testCase.verifyEqual(mark.predict4e(J),hPredict,'AbsTol',testCase.tol);
            testCase.verifyEqual(mark.predict4e(K),mark.gammaP*eta4e(K).^2,'AbsTol',testCase.tol);
        end

        function testPredictorFollowsRefinement(testCase)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            mark = MarkHP;
            mark.predict4e = (1:mesh.nrElems)';
            marked4e = zeros(mesh.nrElems,1);
            marked4e(1:4:end) = 1;
            refined = mesh.refine(marked4e,'NVB');

            mark = mark.initializeMark(refined);

            testCase.verifyEqual(mark.predict4e,refined.parent4e);
            testCase.verifyEqual(mark.marked4e,zeros(refined.nrElems,1));
            testCase.verifyEqual(mark.marked4ed,zeros(refined.nrEdges,1));
        end

        function testAdaptiveLoop(testCase)
            initialMesh = Mesh(Geometry(struct('name','Lshape')));
            mesh = initialMesh;
            mark = MarkBulk;
            for k = 1 : 3
                % indicators concentrated at the reentrant corner (0,0)
                eta4e = 1./(norm4e(mesh.point4e([1 1]/3)) + 0.1);

                mark = mark.mark(struct('mesh',mesh,'degree',1),eta4e);
                refined = mesh.refine(mark.marked4ed);

                % every marked element is split
                nrChildren = accumarray(refined.parent4e,1,[mesh.nrElems,1]);
                testCase.verifyGreaterThanOrEqual(nrChildren(mark.marked4e == 1),2, ...
                    sprintf('step %d',k));
                MeshTest.verifyMesh(testCase,refined,initialMesh,3);
                mesh = refined;
            end
            testCase.verifyGreaterThan(mesh.nrElems,initialMesh.nrElems);
        end

        function testAdaptiveLoop3d(testCase)
            initialMesh = UnitTestCase.cubeMesh('DNDNDN');
            mesh = initialMesh;
            mark = MarkBulk;
            for k = 1 : 3
                eta4e = 1./(norm4e(mesh.point4e([1 1 1]/4)) + 0.1);

                mark = mark.mark(struct('mesh',mesh,'degree',1),eta4e);
                refined = mesh.refine(mark.marked4e);

                testCase.verifyGreaterThan(refined.nrElems,mesh.nrElems,sprintf('step %d',k));
                MeshTest.verifyMesh3d(testCase,refined,initialMesh,1);
                mesh = refined;
            end
        end

    end

    methods (Static)

        function fe = feStruct(degree)
            %feStruct - Mesh and polynomial degree as read by the strategies
            if nargin < 1
                degree = 1;
            end
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = struct('mesh',mesh,'degree',degree, ...
                'degree4e',degree*ones(mesh.nrElems,1));
        end

        function eta = indicators(n)
            %indicators - Distinct positive indicators in scrambled order
            eta = mod(7*(1:n)',n) + 1;
            eta = eta.^1.5/n;
        end

        function marked4ed = edgesOfElements(mesh,marked4e)
            %edgesOfElements - Marks of all edges of the marked elements
            marked4ed = zeros(mesh.nrEdges,1);
            marked4ed(mesh.ed4e(marked4e == 1,:)) = 1;
        end

    end

end

function r = norm4e(x)
%norm4e - Euclidean norm of each row
r = sqrt(sum(x.^2,2));
end
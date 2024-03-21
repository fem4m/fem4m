classdef MeshSurfaceTest < UnitTestCase
    %MeshSurfaceTest - Unit tests of the parametrization of surface meshes
    %   On the unit sphere, the Jacobian of an affine surface mesh is
    %   spanned by the edges of the elements and its derivatives and
    %   curvatures vanish. The isoparametric meshes of higher degree have
    %   Jacobians, second derivatives and normal derivatives that agree with
    %   finite differences of the parametrization and of the unit normal,
    %   the derivative of the unit normal is tangential, the normal
    %   projector is an orthogonal projection onto the tangent plane, the
    %   transformation is a left inverse of the Jacobian, the Pk Hessian
    %   on the curved mesh gives the surface Hessians -n_i*grad(n) and
    %   the surface Laplacian -kappa*n of the coordinates x_i, and the
    %   principal and total curvatures converge to those of the sphere,
    %   1 and 2. The vertices are reprojected after every internal round
    %   of a refine call, so refine('uniform',2) equals two uniform
    %   refinements and graded refinement meets its mesh size.
    %   Coarsening rebuilds the isoparametric map, an empty degree means
    %   degree 1, and surface meshes with symmetry edges survive save and
    %   load.
    %
    %   See also MeshSurface, Pk, problemSurfaceSphere, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degree = struct('P1',1,'P2',2,'P3',3)
        curvedDegree = struct('P2',2,'P3',3)
    end

    properties (Constant)
        ref = [0.2 0.3]  % interior reference point
        h = 1e-5         % step of the finite differences
        fdTol = 1e-7     % tolerance of the finite differences
    end

    methods (Test)

        function testAffineMesh(testCase)
            mesh = MeshSurfaceTest.sphere(1);
            n = mesh.nrElems;
            P = cell(1,3);
            for j = 1 : 3
                P{j} = mesh.c4n(mesh.n4e(:,j),:);
            end

            jacobian = mesh.jacobian(testCase.ref);

            testCase.verifyEqual(reshape(jacobian(:,1,:),n,3),P{2} - P{1});
            testCase.verifyEqual(reshape(jacobian(:,2,:),n,3),P{3} - P{1});
            testCase.verifyEqual(mesh.jacobian4ed(0.5), ...
                mesh.c4n(mesh.n4ed(:,2),:) - mesh.c4n(mesh.n4ed(:,1),:));
            testCase.verifyEqual(mesh.gradJacobian(testCase.ref),zeros(n,3,3));
            testCase.verifyEqual(mesh.DJacobian(testCase.ref),zeros(n,2,3,3));
            testCase.verifyEqual(mesh.surfaceGradNormal4e(testCase.ref),zeros(n,3,3));
            testCase.verifyEqual(mesh.surfaceCurvature4e(testCase.ref),zeros(n,2));
            testCase.verifyEqual(mesh.surfaceTotalCurvature4e(testCase.ref),zeros(n,1));
        end

        function testNormalProjector(testCase,degree)
            mesh = MeshSurfaceTest.sphere(degree);
            n = mesh.nrElems;

            P = mesh.normalProjector(testCase.ref);

            normal = mesh.surfaceNormal4e(testCase.ref);
            testCase.verifySize(P,[n,3,3]);
            testCase.verifyEqual(P,permute(P,[1 3 2]),'AbsTol',testCase.tol,'symmetric');
            testCase.verifyEqual(matMul(P,P),P,'AbsTol',testCase.tol,'projection');
            testCase.verifyEqual(matMul(P,normal),zeros(n,3),'AbsTol',testCase.tol,'P n = 0');
            testCase.verifyEqual(P(:,1,1) + P(:,2,2) + P(:,3,3),2*ones(n,1), ...
                'AbsTol',testCase.tol,'rank 2');
        end

        function testTransformation(testCase,degree)
            mesh = MeshSurfaceTest.sphere(degree);
            n = mesh.nrElems;

            T = mesh.transformation(testCase.ref);

            J = mesh.jacobian(testCase.ref);
            testCase.verifySize(T,[n,3,2]);
            testCase.verifyEqual(matMul(J,T),repmat(reshape(eye(2),1,2,2),n,1,1), ...
                'AbsTol',testCase.tol,'J T = I');
        end

        function testJacobianDerivatives(testCase,curvedDegree)
            mesh = MeshSurfaceTest.sphere(curvedDegree);
            ref = testCase.ref;
            [d1,d2] = MeshSurfaceTest.differences(@(r)mesh.jacobian(r),ref,testCase.h);

            gradJacobian = mesh.gradJacobian(ref);
            DJacobian = mesh.DJacobian(ref);

            % [xx yy xy] derivatives of the parametrization
            tol = testCase.fdTol;
            testCase.verifyEqual(gradJacobian(:,1,:),d1(:,1,:),'AbsTol',tol,'xx');
            testCase.verifyEqual(gradJacobian(:,2,:),d2(:,2,:),'AbsTol',tol,'yy');
            testCase.verifyEqual(gradJacobian(:,3,:),d2(:,1,:),'AbsTol',tol,'xy');
            testCase.verifyEqual(gradJacobian(:,3,:),d1(:,2,:),'AbsTol',tol,'yx');
            % surface gradients of the columns of the Jacobian
            J = mesh.jacobian(ref);
            for i = 1 : 2
                for c = 1 : 3
                    grad = reshape(DJacobian(:,i,c,:),[],3);
                    testCase.verifyEqual(sum(grad.*reshape(J(:,1,:),[],3),2),d1(:,i,c), ...
                        'AbsTol',tol,sprintf('DJacobian(:,%d,%d,:) g1',i,c));
                    testCase.verifyEqual(sum(grad.*reshape(J(:,2,:),[],3),2),d2(:,i,c), ...
                        'AbsTol',tol,sprintf('DJacobian(:,%d,%d,:) g2',i,c));
                end
            end
        end

        function testEdgeJacobian(testCase,curvedDegree)
            mesh = MeshSurfaceTest.sphere(curvedDegree);
            t = 0.3;
            h = testCase.h;

            jacobian4ed = mesh.jacobian4ed(t);

            expected = (mesh.point4ed(t+h) - mesh.point4ed(t-h))/(2*h);
            testCase.verifyEqual(jacobian4ed,expected,'AbsTol',testCase.fdTol);
            testCase.verifyEqual(mesh.tangent4ed(t),expected./sqrt(sum(expected.^2,2)), ...
                'AbsTol',testCase.fdTol);
        end

        function testGradNormal(testCase,curvedDegree)
            mesh = MeshSurfaceTest.sphere(curvedDegree);
            ref = testCase.ref;
            n = mesh.nrElems;

            gradNormal = mesh.surfaceGradNormal4e(ref);

            % the unit normal has tangential derivatives
            normal = mesh.surfaceNormal4e(ref);
            testCase.verifyEqual(reshape(matMul(permute(normal,[1 3 2]),gradNormal),n,3), ...
                zeros(n,3),'AbsTol',testCase.tol,'n'' grad(n) = 0');
            % derivatives along the parametrization
            [d1,d2] = MeshSurfaceTest.differences(@(r)reshape(mesh.surfaceNormal4e(r),n,1,3),ref,testCase.h);
            J = mesh.jacobian(ref);
            testCase.verifyEqual(matMul(gradNormal,reshape(J(:,1,:),n,3)),reshape(d1,n,3), ...
                'AbsTol',testCase.fdTol,'grad(n) g1');
            testCase.verifyEqual(matMul(gradNormal,reshape(J(:,2,:),n,3)),reshape(d2,n,3), ...
                'AbsTol',testCase.fdTol,'grad(n) g2');
        end

        function testRefineRounds(testCase,degree)
            mesh = MeshSurfaceTest.sphere(degree);

            twoRounds = mesh.refine('uniform',2);
            twoCalls = mesh.refine('uniform').refine('uniform');

            testCase.verifyEqual(twoRounds.n4e,twoCalls.n4e);
            testCase.verifyEqual(twoRounds.c4n,twoCalls.c4n,'AbsTol',testCase.tol);
            if degree > 1
                testCase.verifyEqual(twoRounds.pk.c4dof,twoCalls.pk.c4dof,'AbsTol',testCase.tol);
            end
            testCase.verifySize(twoRounds.prolongation,[twoRounds.nrNodes,mesh.nrNodes], ...
                'P1 of the composed rounds');
            H = 0.3;
            graded = mesh.refine('graded',0.5,H,[0 0]);
            length4ed = graded.length4ed;
            testCase.verifyLessThanOrEqual(max(length4ed),H,'graded mesh size');
        end

        function testCoarseSurface(testCase,curvedDegree)
            % coarsening rebuilds the isoparametric map
            p = problemSurfaceSphere();
            mesh = MeshSurfaceTest.sphere(curvedDegree);
            refined = mesh.refine(true(mesh.nrElems,1),'NVB');

            coarse = refined.coarse(true(refined.nrElems,1));

            fresh = MeshSurface(coarse.c4n,coarse.n4e,coarse.n4Db,coarse.n4Nb,coarse.n4Sb, ...
                @p.surfaceFunction,@p.gradSurfaceFunction,curvedDegree);
            testCase.verifyLessThan(coarse.nrElems,refined.nrElems);
            testCase.verifyEqual(coarse.pk.dof4e,fresh.pk.dof4e);
            testCase.verifyEqual(coarse.pk.c4dof,fresh.pk.c4dof,'AbsTol',testCase.tol);
            testCase.verifyEqual(coarse.area4e([1 1]/3),fresh.area4e([1 1]/3),'AbsTol',testCase.tol);
        end

        function testEmptyDegree(testCase)
            p = problemSurfaceSphere();

            mesh = MeshSurface(p.geom,@p.surfaceFunction,@p.gradSurfaceFunction,[]);

            testCase.verifyEqual(mesh.degree,1);
            testCase.verifyEqual(mesh.refine('uniform').nrElems,4*mesh.nrElems);
        end

        function testSymmetrySaveLoad(testCase)
            % the sides of surface meshes are edges: the axis column of a
            % saved mesh is no third node of a face
            p = problemSurfaceSquare();
            g = p.geom;
            n4Sb = [g.n4Db,mod((1:size(g.n4Db,1))',3)];
            mesh = MeshSurface(g.c4n,g.n4e,[],[],n4Sb,@p.surfaceFunction,@p.gradSurfaceFunction);
            file = [tempname,'_mesh'];
            cleanup = onCleanup(@()delete([file,'_*.dat']));

            mesh.save(file);

            testCase.verifyEqual(Geometry(struct('name',file)).n4Sb,n4Sb,'Geometry');
            testCase.verifyEqual(Mesh(file).n4Sb,n4Sb,'Mesh(NAME)');
            testCase.verifyEqual(Mesh(file).c4n,mesh.c4n,'all digits');
            testCase.verifySize(mesh.refine('uniform').n4Sb,[2*size(n4Sb,1),3],'refined');
            % node pairs alone, without the fixed axis
            writematrix(n4Sb(:,1:2),[file,'_n4Sb.dat'],'Delimiter',' ');
            testCase.verifyError(@()Geometry(struct('name',file)),'fem4m:Geometry:invalidN4Sb','Geometry');
            testCase.verifyError(@()MeshSurface(file,@p.surfaceFunction,@p.gradSurfaceFunction), ...
                'fem4m:Mesh:invalidN4Sb','MeshSurface(NAME)');
        end

        function testCurvature(testCase,curvedDegree)
            coarse = MeshSurfaceTest.sphere(curvedDegree).refine('uniform');
            fine = coarse.refine('uniform');

            errors = zeros(1,2);
            meshes = {coarse,fine};
            for k = 1 : 2
                curvature = meshes{k}.surfaceCurvature4e(testCase.ref);
                total = meshes{k}.surfaceTotalCurvature4e(testCase.ref);
                testCase.verifyEqual(total,sum(curvature,2),'AbsTol',testCase.tol, ...
                    'total curvature is the sum of the principal curvatures');
                errors(k) = max(abs([curvature(:) - 1;total/2 - 1]));
            end

            testCase.verifyLessThan(errors(2),errors(1),'convergence');
            testCase.verifyLessThan(errors(2),0.15,'curvature of the unit sphere');
        end

        function testSurfaceHessian(testCase,curvedDegree)
            % the coordinates x_i of the isoparametric surface are Pk
            % functions with surface Hessian -n_i*grad(n) and surface
            % Laplacian -kappa*n_i
            mesh = MeshSurfaceTest.sphere(curvedDegree);
            ref = testCase.ref;
            n = mesh.nrElems;
            fe = Pk(curvedDegree,mesh);
            d2 = fe.getD2Basis([],ref,[]);
            X = reshape(mesh.pk.c4dof(mesh.pk.dof4e,:),n,[],3);
            normal = mesh.surfaceNormal4e(ref);
            gradNormal = mesh.surfaceGradNormal4e(ref);

            % [xx yy zz xy xz yz] as 3 x 3 matrices
            for i = 1 : 3
                H = reshape(sum(d2.*X(:,:,i),2),n,6);
                testCase.verifyEqual(reshape(H(:,[1 4 5 4 2 6 5 6 3]),n,3,3), ...
                    -normal(:,i).*gradNormal,'AbsTol',testCase.tol*max(abs(gradNormal(:))), ...
                    sprintf('Hessian of x_%d',i));
            end
            laplace = reshape(sum(sum(d2(:,:,1:3),3).*X,2),n,3);
            expected = -mesh.surfaceTotalCurvature4e(ref).*normal;
            testCase.verifyEqual(laplace,expected,'AbsTol',testCase.tol*max(abs(expected(:))), ...
                'Laplacian');
        end

    end

    methods (Static)

        function mesh = sphere(degree)
            %sphere - Isoparametric surface mesh of the unit sphere
            p = problemSurfaceSphere();
            mesh = MeshSurface(p.geom,@p.surfaceFunction,@p.gradSurfaceFunction,degree);
        end

        function [d1,d2] = differences(f,ref,h)
            %differences - Central differences of F in the reference point
            d1 = (f(ref + [h 0]) - f(ref - [h 0]))/(2*h);
            d2 = (f(ref + [0 h]) - f(ref - [0 h]))/(2*h);
        end

    end

end
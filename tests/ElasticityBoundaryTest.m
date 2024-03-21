classdef ElasticityBoundaryTest < UnitTestCase
    %ElasticityBoundaryTest - Unit tests of the elasticity boundary conditions
    %   The traction data g(x,n) = sigma(x)*n of a stress field sigma in
    %   the Arnold-Winther type elements AW1, AWNC1 and HZ1 gives, by the
    %   scaled boundary moments and corner values of interpolateBoundary,
    %   the degrees of freedom of the projection of sigma on the Neumann
    %   boundary and zero elsewhere. On a symmetry boundary the elements
    %   fix the tangential traction moments: they vanish in the
    %   projection of a stress field with sigma_12 = 0 and not for
    %   sigma_12 ~= 0. AWH adds three nodal dofs at each of the four
    %   non-parallel Neumann corners of the square, and the traction data
    %   still gives the boundary dofs of the projection. AWk(3) is the
    %   space of AW1 and integrates the boundary moments at least with
    %   the quadrature of AW1, so non-polynomial traction data give the
    %   same boundary function in both elements.
    %
    %   On the faces of a 3d symmetry boundary PkVector fixes the
    %   displacement component of the axis, and the displacement
    %   u = (x, 2y, 3z) gliding on the faces x = 0 and z = 0 solves the
    %   discrete problem exactly.
    %
    %   See also AW1, AWNC1, HZ1, AWH, PkVector, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        element = struct('AW1',@(m)AW1(3,m),'AWNC1',@(m)AWNC1(2,m),'HZ1',@(m)HZ1(3,m))
    end

    methods (Test)

        function testTractionData(testCase,element)
            mesh = PatchTestCase.patchMesh(struct('name','Square','bdry','N'));
            fe = element(mesh);
            sigma = @(z)[z(:,1) + 1, z(:,1) - z(:,2), z(:,1) - z(:,2), z(:,2) + 2];
            g = @(z,n)ElasticityBoundaryTest.traction(sigma(z),n);
            x = fe.project(sigma);

            val = fe.interpolateBoundary(g);

            set = val ~= 0;
            set(fe.fixedNodes) = true;
            testCase.assertGreaterThan(nnz(val),0,'boundary data');
            testCase.verifyEqual(val(set),x(set),'AbsTol',testCase.tol,'boundary degrees of freedom');
        end

        function testSymmetryBoundary(testCase,element)
            mesh = PatchTestCase.patchMesh(struct('name','Square','bdry','DSDS'));
            fe = element(mesh);
            % without Neumann boundary only the symmetry dofs are fixed
            fixed = fe.fixedNodes(:);
            onSb = unique(fe.dof4ed(mesh.ed4Sb,:));
            testCase.assertNotEmpty(fixed,'symmetry degrees of freedom');
            testCase.verifyEmpty(setdiff(fixed,onSb),'fixed dofs on the symmetry boundary');

            x = fe.project(@(z)[z(:,1) + 1, 0*z(:,1), 0*z(:,1), z(:,2) + 2]);
            testCase.verifyEqual(x(fixed),zeros(size(fixed)),'AbsTol',testCase.tol,'sigma_12 = 0');
            testCase.verifyGreaterThan(max(abs(x(setdiff(onSb,fixed)))),0.1,'normal traction');

            x = fe.project(@(z)[z(:,1) + 1, 1 + z(:,2), 1 + z(:,2), z(:,2) + 2]);
            testCase.verifyGreaterThan(max(abs(x(fixed))),0.1,'sigma_12 ~= 0');
        end

        function testAWkBoundaryDataOfAW1(testCase)
            mesh = PatchTestCase.patchMesh(struct('name','Square','bdry','N'));
            % oscillating data, which a Gauss rule with fewer points than
            % that of AW1 integrates with an error far above round-off
            sigma = @(z)[sin(6*z(:,1)), cos(6*z(:,2)), cos(6*z(:,2)), exp(z(:,1).*z(:,2))];
            g = @(z,n)ElasticityBoundaryTest.traction(sigma(z),n);
            aw1 = AW1(3,mesh);
            awk = AWk(3,mesh);
            ref = [0.2 0.3];

            val = awk.getValue([],ref,[],awk.interpolateBoundary(g));

            expected = aw1.getValue([],ref,[],aw1.interpolateBoundary(g));
            testCase.verifyEqual(val,expected,'AbsTol',1e-12*max(abs(expected(:))));
        end

        function testHybridCorners(testCase)
            % every corner of the square has two triangles after NVB
            mesh = Mesh(struct('name','Square','bdry','N')).refine('uniform','nvb');
            fe = AWH(3,mesh);
            sigma = @(z)[z(:,1) + 1, z(:,1) - z(:,2), z(:,1) - z(:,2), z(:,2) + 2];
            g = @(z,n)ElasticityBoundaryTest.traction(sigma(z),n);
            x = fe.project(sigma);

            val = fe.interpolateBoundary(g);

            testCase.verifyEqual(fe.hybridNodes,4,'hybrid corners');
            testCase.verifyEqual(fe.nrDoF,AWk(3,mesh).nrDoF + 12,'hybrid dofs');
            set = val ~= 0;
            set(fe.fixedNodes) = true;
            testCase.verifyEqual(val(set),x(set),'AbsTol',testCase.tol,'boundary degrees of freedom');
        end

    end

    methods (Test)

        function testGlidingBoundary3d(testCase)
            mesh = Mesh3d(Geometry(struct('name','Cube','bdry','DDDSDS')));
            mesh = mesh.refine('uniform');
            fe = PkVector(2,mesh);
            u = @(z)[z(:,1), 2*z(:,2), 3*z(:,3)];
            % the faces x = 0 (axis 1) and z = 0 (axis 3)
            for k = [1 3]
                onSb = (k-1)*fe.fe{1}.nrDoF + unique(fe.fe{k}.dof4f(mesh.f4Sb(mesh.n4Sb(:,4)==k),:));
                testCase.verifyTrue(all(ismember(onSb,fe.fixedNodes)),'fixed component');
            end
            onSb2 = fe.fe{1}.nrDoF + unique(fe.fe{2}.dof4f(mesh.f4Sb,:));
            onDb2 = fe.fe{1}.nrDoF + unique(fe.fe{2}.dof4f(mesh.f4Db,:));
            testCase.verifyEmpty(intersect(fe.fixedNodes,setdiff(onSb2,onDb2)), ...
                'free tangential component');

            A = BilinearForm('eps(u):eps(v) dx').assemble(fe);
            x = zeros(fe.nrDoF,1);
            uInterpolated = fe.interpolate(u);
            x(fe.fixedNodes) = uInterpolated(fe.fixedNodes);
            b = -A*x;
            x(fe.freeNodes) = A(fe.freeNodes,fe.freeNodes)\b(fe.freeNodes);

            err = ErrorForm('(exact(x)-u)^2 dx').assemble(x,fe,'exact',u);
            testCase.verifyLessThan(err,testCase.tol,'gliding solution');
        end

    end

    methods (Static)

        function val = traction(sigma,n)
            %traction - sigma*n for the matrices sigma stored column by column
            val = [sigma(:,1).*n(:,1) + sigma(:,3).*n(:,2), ...
                sigma(:,2).*n(:,1) + sigma(:,4).*n(:,2)];
        end

    end

end
classdef FormTest < UnitTestCase
    %FormTest - Unit tests of the assembly of bilinear and linear forms
    %   The stiffness matrix 'grad(u)*grad(v) dx' without boundary
    %   conditions is symmetric positive semidefinite with the constants
    %   in its kernel, and the mass matrix 'u*v dx' is symmetric positive
    %   definite with one'*M*one equal to the area of the domain for the
    %   coefficients one of the constant function. For Lagrange elements
    %   one is the vector of ones, so that the sum of all entries of M is
    %   the area (volume) of the domain. For the
    %   interpolant x of a polynomial f in the discrete space, x'*A*x and
    %   x'*M*x are the exact integrals of |grad(f)|^2 and f^2, and the
    %   load vector of 'f(x) v dx' equals M*x. The boundary mass matrix
    %   'u*v dsd' sums up to the length of the Dirichlet boundary, and the
    %   error form '(u)^2' of the constant one on the boundary domains is
    %   the square root of their length.
    %
    %   Assembled on a subdomain of every other element, the mass matrix
    %   and the load vector of the constant one sum up to the area of the
    %   subdomain, and the forms on the sides of the Dirichlet boundary
    %   as subdomain of 'ds' equal those on 'dsd'. For the elementwise
    %   constant function g = 1 + T/nrElems in the discontinuous spaces,
    %   x'*b for the load vectors b of 'f(x)*[v] ds' and 'f(x)*{v} dsi'
    %   with f = 1 + x_1 are the integrals of f times the jump (the trace
    %   on the boundary) and the average of g.
    %
    %   The parser does not replace the letters of operator keywords and
    %   coefficient names by the variables, e.g. for the variables r, c,
    %   l, o, t in curl and rot, v in div and dev, or u in mu(x), also
    %   without parentheses as in u*divv or curlu*curl v. The
    %   jump of the first operator of a variable applies to all of its
    %   operators, both for [a]-[b] and [a-b]. The curl-type operators
    %   satisfy |curl u| = |grad u| and
    %   curl(curl(u)):curl(curl(v)) = D^2 u:D^2 v for scalar u, v in 2d,
    %   and (curl u,v) = (u,rot v) in 2d and (curl u,v) = (u,curl v) in
    %   3d for u and v vanishing on the boundary.
    %
    %   All argument variants of BilinearForm.assemble (degree, subdomain
    %   and named coefficients after FE1 or FE2) give the matrix of the
    %   fully specified call, and the local element and side matrices sum
    %   up to the matrix, also for negative terms. The forms on the
    %   symmetry boundary dss (in 2d and 3d) and on the faces of a 3d
    %   mesh sum up to the length and area of these sides. The estimator of the difference
    %   of two functions u - v, of a weighted term, and of symmetric
    %   matrix fields with a matrix coefficient or a source function in
    %   Voigt notation give the corresponding element norms.
    %
    %   Form.sideOf gives the edge data of 2d, surface and polygonal
    %   meshes and the face data of Mesh3d under one name, and the side
    %   dofs and degrees of the elements, [] without them, also of the
    %   edge and face spaces.
    %
    %   See also BilinearForm, LinearForm, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        element = struct( ...
            'P1',struct('name','Pk','degree',1), ...
            'P2',struct('name','Pk','degree',2), ...
            'P3',struct('name','Pk','degree',3), ...
            'CR1',struct('name','CRk','degree',1), ...
            'DGP2',struct('name','DGPk','degree',2), ...
            'DGD2',struct('name','DGDk','degree',2))
        lagrange = struct('P2',2,'P3',3,'P4',4)
        dgElement = struct( ...
            'DGP1',struct('name','DGPk','degree',1), ...
            'DGP2',struct('name','DGPk','degree',2), ...
            'DGD0',struct('name','DGDk','degree',0), ...
            'DGD2',struct('name','DGDk','degree',2))
        element3d = struct('P1',1,'P2',2,'P3',3)
        flux = struct( ...
            'RT1',struct('name','RTk','degree',1), ...
            'RT2',struct('name','RTk','degree',2), ...
            'BDM1',struct('name','BDMk','degree',1))
        parserCase = struct( ...
            'curlR',struct('str','curl(r)*curl(s) dx','varStr',{{'r','s'}},'ops',{{'curlw_T','curlw_T'}}), ...
            'curlC',struct('str','curl(c)*curl(q) dx','varStr',{{'c','q'}},'ops',{{'curlw_T','curlw_T'}}), ...
            'curlL',struct('str','curl(l)*curl(q) dx','varStr',{{'l','q'}},'ops',{{'curlw_T','curlw_T'}}), ...
            'rotO',struct('str','rot(o)*rot(q) dx','varStr',{{'o','q'}},'ops',{{'rotw_T','rotw_T'}}), ...
            'rotT',struct('str','rot(t)*rot(q) dx','varStr',{{'t','q'}},'ops',{{'rotw_T','rotw_T'}}), ...
            'divV',struct('str','div(v)*div(u) dx','varStr',{{'v','u'}},'ops',{{'divw_T','divw_T'}}), ...
            'divNoParentheses',struct('str','u*divv dx','varStr',{{}},'ops',{{'w_T','divw_T'}}), ...
            'curlNoParentheses',struct('str','curlu*curl v dx','varStr',{{}},'ops',{{'curlw_T','curlw_T'}}), ...
            'graddivNoParentheses',struct('str','graddivsigma*divtau dx','varStr',{{'sigma','tau'}},'ops',{{'graddivw_T','divw_T'}}), ...
            'dev',struct('str','dev(u):dev(v) dx','varStr',{{}},'ops',{{'devw_T','devw_T'}}), ...
            'coefficient',struct('str','mu(x)*u*v dx','varStr',{{}},'ops',{{'w_T','w_T'}}), ...
            'curlCurl',struct('str','curl(curl(u)):v dx','varStr',{{}},'ops',{{'curlcurlw_T','w_T'}}), ...
            'curlRot',struct('str','curl(rot(u))*curlrot(v) dx','varStr',{{}},'ops',{{'curlrotw_T','curlrotw_T'}}), ...
            'rotRot',struct('str','Rot(Rot(sigma))*rotrot(tau) dx','varStr',{{'sigma','tau'}},'ops',{{'rotrotw_T','rotrotw_T'}}), ...
            'side',struct('str','u*v ds','varStr',{{}},'ops',{{'w_S','w_S'}}), ...
            'sideJump',struct('str','[u]*[v] dsi','varStr',{{}},'ops',{{'jumpw_S','jumpw_S'}}), ...
            'sideJumps',struct('str','([grad(u)*n]-[u])*[v] dsi','varStr',{{}},'ops',{{'jumpgradwn_S','jumpw_S','jumpw_S'}}), ...
            'sideJumpOfDifference',struct('str','[grad(u)*n-u]*[v] dsi','varStr',{{}},'ops',{{'jumpgradwn_S','jumpw_S','jumpw_S'}}), ...
            'sideAverage',struct('str','{grad(u)*n}*v dsb','varStr',{{}},'ops',{{'averagegradwn_S','w_S'}}), ...
            'sym',struct('str','sym(u):sym(v) dx','varStr',{{}},'ops',{{'symw_T','symw_T'}}), ...
            'symCurl',struct('str','sym(curl(u)):sym(curl(v)) dx','varStr',{{}},'ops',{{'symcurlw_T','symcurlw_T'}}), ...
            'divDiv',struct('str','divdiv(u)*div(div(v)) dx','varStr',{{}},'ops',{{'divdivw_T','divdivw_T'}}), ...
            'divD2',struct('str','div(d2(u))*grad(v) dx','varStr',{{}},'ops',{{'divd2w_T','gradw_T'}}), ...
            'sideD2DivNormal',struct('str','d2(u)*n*div(v)*n ds','varStr',{{}},'ops',{{'d2wn_S','divwn_S'}}), ...
            'sideDivD2Normal',struct('str','div(d2(u))*n*v ds','varStr',{{}},'ops',{{'divd2wn_S','w_S'}}), ...
            'sideJumpSymCurl',struct('str','[sym(curl(u))]*t*[v]*t dsi','varStr',{{}},'ops',{{'jumpsymcurlwt_S','jumpwt_S'}}), ...
            'sideJumpD2',struct('str','[d2(u)]*t*[div(d2(v))]*t dsi','varStr',{{}},'ops',{{'jumpd2wt_S','jumpdivd2wt_S'}}), ...
            'sideJumpDiv',struct('str','[div(u)]*n*[div(v)]*t dsi','varStr',{{}},'ops',{{'jumpdivwn_S','jumpdivwt_S'}}), ...
            'sideAverageD2',struct('str','{div(d2(u))}*t*{d2(v)}*t dsi','varStr',{{}},'ops',{{'averagedivd2wt_S','averaged2wt_S'}}), ...
            'sideAverageDiv',struct('str','{div(u)}*n*{div(v)}*t dsi','varStr',{{}},'ops',{{'averagedivwn_S','averagedivwt_S'}}))
    end

    methods (Test)

        function testStiffnessMatrix(testCase,element)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = feval(element.name,element.degree,mesh);

            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);

            FormTest.verifySymmetric(testCase,A);
            FormTest.verifySemidefinite(testCase,A);
            one = fe.project(@(x)ones(size(x,1),1));
            testCase.verifyEqual(A*one,zeros(fe.nrDoF,1), ...
                'AbsTol',testCase.tol,'constants in the kernel');
        end

        function testMassMatrix(testCase,element)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = feval(element.name,element.degree,mesh);

            M = BilinearForm('u * v dx').assemble(fe);

            FormTest.verifySymmetric(testCase,M);
            FormTest.verifyDefinite(testCase,M);
            one = fe.project(@(x)ones(size(x,1),1));
            testCase.verifyEqual(one'*M*one,3,'AbsTol',testCase.tol,'area');
        end

        function testStiffnessMatrix3d(testCase,element3d)
            fe = Pk3d(element3d,UnitTestCase.cubeMesh());

            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);

            FormTest.verifySymmetric(testCase,A);
            FormTest.verifySemidefinite(testCase,A);
            testCase.verifyEqual(A*ones(fe.nrDoF,1),zeros(fe.nrDoF,1), ...
                'AbsTol',testCase.tol,'constants in the kernel');
        end

        function testMassMatrix3d(testCase,element3d)
            fe = Pk3d(element3d,UnitTestCase.cubeMesh());

            M = BilinearForm('u * v dx').assemble(fe);

            FormTest.verifySymmetric(testCase,M);
            FormTest.verifyDefinite(testCase,M);
            testCase.verifyEqual(full(sum(M(:))),1,'AbsTol',testCase.tol,'volume');
        end

        function testFluxMassMatrix(testCase,flux)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = feval(flux.name,flux.degree,mesh);

            M = BilinearForm('sigma * tau dx',{'sigma','tau'}).assemble(fe);

            FormTest.verifySymmetric(testCase,M);
            FormTest.verifyDefinite(testCase,M);
        end

        function testEnergyOfInterpolant(testCase,lagrange)
            % f = x^2 - y on the unit square: int |grad f|^2 = 4/3 + 1 and
            % int f^2 = 1/5 - 1/3 + 1/3 = 1/5
            mesh = PatchTestCase.patchMesh(struct('name','Square'));
            fe = Pk(lagrange,mesh);
            x = fe.c4dof(:,1).^2 - fe.c4dof(:,2);

            A = BilinearForm('grad(u) * grad(v) dx').assemble(fe);
            M = BilinearForm('u * v dx').assemble(fe);

            testCase.verifyEqual(x'*A*x,7/3,'AbsTol',testCase.tol,'H1 seminorm');
            testCase.verifyEqual(x'*M*x,1/5,'AbsTol',testCase.tol,'L2 norm');
        end

        function testLoadVector(testCase,lagrange)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = Pk(lagrange,mesh);
            f = UnitTestCase.polynomial(lagrange);

            M = BilinearForm('u * v dx').assemble(fe);
            b = LinearForm('f(x) v dx').assemble(fe,'f',f);

            testCase.verifyEqual(b,M*f(fe.c4dof),'AbsTol',testCase.tol);
        end

        function testBoundaryMassMatrix(testCase,lagrange)
            mesh = PatchTestCase.patchMesh(struct('name','Square','bdry','DNDN'));
            fe = Pk(lagrange,mesh);

            MD = BilinearForm('u * v dsd').assemble(fe);
            MN = BilinearForm('u * v dsn').assemble(fe);

            FormTest.verifySymmetric(testCase,MD);
            FormTest.verifySemidefinite(testCase,MD);
            testCase.verifyEqual(full(sum(MD(:))),2,'AbsTol',testCase.tol,'Dirichlet length');
            testCase.verifyEqual(full(sum(MN(:))),2,'AbsTol',testCase.tol,'Neumann length');
            testCase.verifyLessThan(full(max(abs(MD(fe.freeNodes,:)),[],'all')), ...
                testCase.tol,'only Dirichlet dofs');
        end

        function testBoundaryErrorForm(testCase,lagrange)
            mesh = PatchTestCase.patchMesh(struct('name','Square','bdry','DNDN'));
            fe = Pk(lagrange,mesh);
            one = ones(fe.nrDoF,1);

            testCase.verifyEqual(ErrorForm('(u)^2 dsb').assemble(one,fe),2, ...
                'AbsTol',testCase.tol,'whole boundary');
            testCase.verifyEqual(ErrorForm('(u)^2 dsd').assemble(one,fe),sqrt(2), ...
                'AbsTol',testCase.tol,'Dirichlet boundary');
            testCase.verifyEqual(ErrorForm('(u)^2 dsn').assemble(one,fe),sqrt(2), ...
                'AbsTol',testCase.tol,'Neumann boundary');
        end

        function testAssembleArguments(testCase)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = Pk(2,mesh);
            sub = 1:2:mesh.nrElems;
            a = @(x)1 + x(:,1);
            m = BilinearForm('u * v dx');
            form = BilinearForm('alpha(x) u * v dx');
            M = m.assemble(fe,fe,5,sub);
            A = form.assemble(fe,fe,5,sub,'alpha',a);
            Aall = form.assemble(fe,fe,5,[],'alpha',a);

            testCase.verifyEqual(m.assemble(fe,5,sub),M,'AbsTol',testCase.tol,'(fe,degree,sub)');
            testCase.verifyEqual(m.assemble(fe,sub),m.assemble(fe,fe,[],sub), ...
                'AbsTol',testCase.tol,'(fe,sub)');
            testCase.verifyEqual(m.assemble(fe,5),m.assemble(fe,fe,5), ...
                'AbsTol',testCase.tol,'(fe,degree)');
            testCase.verifyEqual(form.assemble(fe,'alpha',a),form.assemble(fe,fe,[],[],'alpha',a), ...
                'AbsTol',testCase.tol,'(fe,name,value)');
            testCase.verifyEqual(form.assemble(fe,5,'alpha',a),Aall, ...
                'AbsTol',testCase.tol,'(fe,degree,name,value)');
            testCase.verifyEqual(form.assemble(fe,sub,'alpha',a),form.assemble(fe,fe,[],sub,'alpha',a), ...
                'AbsTol',testCase.tol,'(fe,sub,name,value)');
            testCase.verifyEqual(form.assemble(fe,fe,'alpha',a),form.assemble(fe,fe,[],[],'alpha',a), ...
                'AbsTol',testCase.tol,'(fe,fe,name,value)');
            testCase.verifyEqual(form.assemble(fe,fe,5,'alpha',a),Aall, ...
                'AbsTol',testCase.tol,'(fe,fe,degree,name,value)');
            testCase.verifyEqual(form.assemble(fe,5,sub,'alpha',a),A, ...
                'AbsTol',testCase.tol,'(fe,degree,sub,name,value)');
        end

        function testLocalMatrices(testCase)
            mesh = PatchTestCase.patchMesh(struct('name','Square','bdry','DNDN'));
            fe = Pk(2,mesh);

            [A,localA] = BilinearForm('grad(u) * grad(v) dx - u * v dx').assemble(fe);
            [MD,localMD] = BilinearForm('u * v dsd').assemble(fe);

            [I,J] = localDoFtoGlobalDoF(fe.dof4e);
            testCase.verifyEqual(sparse(I(:),J(:),localA(:),fe.nrDoF,fe.nrDoF),A, ...
                'AbsTol',testCase.tol,'element matrices');
            % side matrices act on the dofs of the first neighbor
            [I,J] = localDoFtoGlobalDoF(fe.dof4e(mesh.e4ed(mesh.ed4Db,1),:));
            local = localMD(mesh.ed4Db,:,:);
            testCase.verifyEqual(sparse(I(:),J(:),local(:),fe.nrDoF,fe.nrDoF),MD, ...
                'AbsTol',testCase.tol,'side matrices');
            others = setdiff(1:mesh.nrEdges,mesh.ed4Db);
            testCase.verifyEqual(nnz(localMD(others,:,:)),0,'no other sides');
        end

        function testSymmetryBoundary(testCase)
            mesh = PatchTestCase.patchMesh(struct('name','Square','bdry','DSNS'));
            fe = Pk(2,mesh);
            one = ones(fe.nrDoF,1);
            length4ed = mesh.length4ed;
            lengthS = sum(length4ed(mesh.ed4Sb));
            testCase.assertGreaterThan(lengthS,0,'symmetry boundary');

            M = BilinearForm('u * v dss').assemble(fe);
            eta4e = EstimatorForm('([u])^2 dss').assemble(one,fe);

            testCase.verifyEqual(one'*M*one,lengthS,'AbsTol',testCase.tol,'mass matrix');
            testCase.verifyEqual(sum(eta4e.^2),lengthS,'AbsTol',testCase.tol,'estimator');
            testCase.verifyEqual(full(sum(BilinearForm('u * v dsb').assemble(fe),'all')),4, ...
                'AbsTol',testCase.tol,'whole boundary');
        end

        function testSymmetryBoundary3d(testCase)
            mesh = Mesh3d(Geometry(struct('name','Cube','bdry','DSNSDS')));
            mesh = mesh.refine('uniform');
            fe = Pk3d(2,mesh);
            one = ones(fe.nrDoF,1);
            area4f = mesh.area4f;
            testCase.assertEqual(sum(area4f(mesh.f4Sb)),3,'AbsTol',testCase.tol,'symmetry boundary');

            M = BilinearForm('u * v dss').assemble(fe);
            b = LinearForm('f(x) * v dss').assemble(fe,'f',@(x)ones(size(x,1),1));
            eta4e = EstimatorForm('(u)^2 dss').assemble(one,fe);
            MInterior = BilinearForm('u * v dsi').assemble(fe);

            testCase.verifyEqual(one'*M*one,3,'AbsTol',testCase.tol,'mass matrix');
            testCase.verifyEqual(sum(b),3,'AbsTol',testCase.tol,'load vector');
            testCase.verifyEqual(sum(eta4e.^2),3,'AbsTol',testCase.tol,'estimator');
            testCase.verifyEqual(ErrorForm('(u)^2 dss').assemble(one,fe),sqrt(3), ...
                'AbsTol',testCase.tol,'error form');
            testCase.verifyEqual(one'*MInterior*one,sum(area4f)-6,'AbsTol',testCase.tol, ...
                'interior faces');
            testCase.verifyEqual(full(sum(BilinearForm('u * v dsb').assemble(fe),'all')),6, ...
                'AbsTol',testCase.tol,'whole boundary');
        end

        function testSideForms3d(testCase)
            mesh = UnitTestCase.cubeMesh();
            fe = Pk3d(2,mesh);
            one = ones(fe.nrDoF,1);

            M = BilinearForm('u * v dsd').assemble(fe);
            etaD = EstimatorForm('(u)^2 dsd').assemble(one,fe);
            etaS = EstimatorForm('(u)^2 ds').assemble(one,fe);
            etaJump = EstimatorForm('([u])^2 ds').assemble(one,fe);

            testCase.verifyEqual(one'*M*one,6,'AbsTol',testCase.tol,'Dirichlet boundary');
            testCase.verifyEqual(sum(etaD.^2),6,'AbsTol',testCase.tol,'estimator on the boundary');
            % interior faces count half to both neighbors
            testCase.verifyEqual(sum(etaS.^2),sum(mesh.area4f),'AbsTol',testCase.tol, ...
                'estimator on all faces');
            % the jump is the trace on the boundary and vanishes inside
            testCase.verifyEqual(sum(etaJump.^2),6,'AbsTol',testCase.tol,'jump estimator');
            testCase.verifyEqual(ErrorForm('(u)^2 dsd').assemble(one,fe),sqrt(6), ...
                'AbsTol',testCase.tol,'error form on the boundary');
        end

        function testSideOf(testCase)
            % the sides are the edges of 2d, surface and polygonal meshes
            % and the faces of Mesh3d, decided by the class and not by the
            % number of coordinates (a MeshSurface has three)
            names = {'n4s','e4s','s4e','nrSides','s4Db','s4Nb','s4Sb','normal4s'};
            names2d = {'n4ed','e4ed','ed4e','nrEdges','ed4Db','ed4Nb','ed4Sb','normal4ed'};
            names3d = {'n4f','e4f','f4e','nrFaces','f4Db','f4Nb','f4Sb','normal4f'};
            p = problemSurfaceSquare();
            meshes = {Mesh(Geometry(struct('name','Square','bdry','SDSN'))).refine('uniform'), ...
                MeshSurface(p.geom,@p.surfaceFunction,@p.gradSurfaceFunction,[]), ...
                MeshPolyH(Geometry(struct('name','Square','type','H','bdry','SDSN'))), ...
                UnitTestCase.cubeMesh('DSNSDS')};
            for k = 1 : numel(meshes)
                mesh = meshes{k};
                length4ed = mesh.length4ed;
                if isa(mesh,'Mesh3d')
                    side = 'Face';
                    native = names3d;
                    diam4s = max(length4ed(mesh.ed4f),[],2);
                else
                    side = 'Edge';
                    native = names2d;
                    diam4s = length4ed;
                end
                testCase.verifyEqual(Form.sideOf(mesh),side,class(mesh));
                values = cell(1,numel(names));
                [values{:}] = Form.sideOf(mesh,names{:});
                for j = 1 : numel(names)
                    testCase.verifyEqual(values{j},mesh.(native{j}),[class(mesh),' ',names{j}]);
                end
                testCase.verifyEqual(Form.sideOf(mesh,'diam4s'),diam4s,[class(mesh),' diam4s']);
            end
            % side dofs and degrees of the elements, [] without them, also
            % of the edge and face spaces DGDkEdge and DGDkFace, which are
            % no FE subclasses
            mesh = meshes{1};
            fe = HPk(3,mesh);
            [dof4s,degree4s,e4s] = Form.sideOf(fe,'dof4s','degree4s','e4s');
            testCase.verifyEqual(Form.sideOf(fe),'Edge','HPk side');
            testCase.verifyEqual(dof4s,fe.dof4ed,'HPk dof4s');
            testCase.verifyEqual(degree4s,fe.degree4ed,'HPk degree4s');
            testCase.verifyEqual(e4s,mesh.e4ed,'mesh data of FE.mesh');
            fe = Pk(2,mesh);
            testCase.verifyEqual(Form.sideOf(fe,'dof4s'),fe.dof4ed,'Pk dof4s');
            testCase.verifyEmpty(Form.sideOf(fe,'degree4s'),'Pk degree4s');
            testCase.verifyEmpty(Form.sideOf(DGDk(1,mesh),'dof4s'),'DGDk dof4s');
            fe = DGDkEdge(1,mesh);
            testCase.verifyEqual(Form.sideOf(fe,'dof4s'),fe.dof4ed,'DGDkEdge dof4s');
            testCase.verifyEmpty(Form.sideOf(fe,'degree4s'),'DGDkEdge degree4s');
            mesh = meshes{4};
            fe = Pk3d(3,mesh);
            testCase.verifyEqual(Form.sideOf(fe),'Face','Pk3d side');
            testCase.verifyEqual(Form.sideOf(fe,'dof4s'),fe.dof4f,'Pk3d dof4s');
            fe = DGDkFace(1,mesh);
            testCase.verifyEqual(Form.sideOf(fe,'dof4s'),fe.dof4f,'DGDkFace dof4s');
            testCase.verifyEmpty(Form.sideOf(fe,'degree4s'),'DGDkFace degree4s');
            fe = DGPkFace(1,mesh);
            testCase.verifyEqual(Form.sideOf(fe,'dof4s'),fe.dof4f,'DGPkFace dof4s');
            testCase.verifyEmpty(Form.sideOf(fe,'degree4s'),'DGPkFace degree4s');
        end

        function testEstimatorTwoSpaces(testCase)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = Pk(2,mesh);
            c = fe.project(UnitTestCase.polynomial(2));
            area4e = mesh.area4e;

            eta = EstimatorForm('(u - v)^2 dx').assemble(c,fe,c + 1,fe);
            etaDegree = EstimatorForm('(u - v)^2 dx').assemble(c,fe,c + 1,fe,4);

            testCase.verifyEqual(eta,sqrt(area4e),'AbsTol',testCase.tol,'constant difference');
            testCase.verifyEqual(etaDegree,sqrt(area4e),'AbsTol',testCase.tol,'with degree');
            testCase.verifyWarning(@()EstimatorForm('(u)^2 dx').assemble(c,fe,c + 1,fe), ...
                'fem4m:EstimatorForm:defaultVarStr');
        end

        function testEstimatorWeight(testCase)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = Pk(1,mesh);
            one = ones(fe.nrDoF,1);

            eta = EstimatorForm('(u)^2 w(x) dx').assemble(one,fe,'w',@(x)2 + 0*x(:,1));

            testCase.verifyEqual(eta,sqrt(2*mesh.area4e),'AbsTol',testCase.tol);
        end

        function testEstimatorTensorField(testCase)
            % constant symmetric matrix fields S with the coefficient
            % 2*I in Voigt notation and as source function
            S2 = [1 2;2 3];
            S3 = [1 2 3;2 4 5;3 5 6];
            for S = {S2,S3}
                dim = size(S{1},1);
                if dim == 2
                    mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
                    voigt = S{1}([1 4 3]);
                else
                    mesh = UnitTestCase.cubeMesh();
                    voigt = S{1}([1 5 9 4 7 8]);
                end
                fe = DGPkSymMatrix(1,mesh);
                c = fe.project(@(z)ones(size(z,1),1)*S{1}(:)');
                nrVoigt = numel(voigt);

                eta = EstimatorForm('(alpha(x) u)^2 dx').assemble(c,fe,'alpha',2*eye(nrVoigt));
                etaSource = EstimatorForm('(f(x) - u)^2 dx').assemble(c,fe, ...
                    'f',@(z)ones(size(z,1),1)*voigt);

                testCase.verifyEqual(eta,2*norm(S{1},'fro')*sqrt(mesh.area4e), ...
                    'AbsTol',testCase.tol,sprintf('coefficient in %dd',dim));
                testCase.verifyEqual(etaSource,zeros(mesh.nrElems,1),'AbsTol',testCase.tol, ...
                    sprintf('source function in %dd',dim));
            end
        end

        function testSubdomainMassMatrix(testCase,element)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = feval(element.name,element.degree,mesh);
            subdomain = 1:2:mesh.nrElems;
            area4e = mesh.area4e;

            M = BilinearForm('u * v dx').assemble(fe,fe,[],subdomain);

            one = fe.project(@(x)ones(size(x,1),1));
            testCase.verifyEqual(one'*M*one,sum(area4e(subdomain)), ...
                'AbsTol',testCase.tol,'area of the subdomain');
        end

        function testSubdomainLoadVector(testCase,element)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = feval(element.name,element.degree,mesh);
            subdomain = 1:2:mesh.nrElems;
            area4e = mesh.area4e;

            b = LinearForm('f(x) v dx').assemble(fe,2*fe.degree,subdomain, ...
                'f',@(x)ones(size(x,1),1));

            one = fe.project(@(x)ones(size(x,1),1));
            testCase.verifyEqual(one'*b,sum(area4e(subdomain)), ...
                'AbsTol',testCase.tol,'area of the subdomain');
        end

        function testBoundarySubdomain(testCase,lagrange)
            mesh = PatchTestCase.patchMesh(struct('name','Square','bdry','DNDN'));
            fe = Pk(lagrange,mesh);
            g = @(x)1 + x(:,1);

            M = BilinearForm('u * v ds').assemble(fe,fe,[],mesh.ed4Db);
            b = LinearForm('g(x) v ds').assemble(fe,2*lagrange,mesh.ed4Db,'g',g);

            testCase.verifyEqual(M,BilinearForm('u * v dsd').assemble(fe), ...
                'AbsTol',testCase.tol,'mass matrix');
            testCase.verifyEqual(b,LinearForm('g(x) v dsd').assemble(fe,'g',g), ...
                'AbsTol',testCase.tol,'load vector');
        end

        function testJumpLoadVector(testCase,dgElement)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = feval(dgElement.name,dgElement.degree,mesh);
            [x,g4e] = FormTest.elementwiseConstant(fe);
            [weight,g1,g2] = FormTest.sideIntegrals(mesh,g4e);

            b = LinearForm('f(x)*[v] ds').assemble(fe,'f',@(x)1 + x(:,1));

            testCase.verifyEqual(x'*b,sum(weight.*(g1 - g2)), ...
                'AbsTol',testCase.tol,'integral of f times the jump');
        end

        function testAverageLoadVector(testCase,dgElement)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            fe = feval(dgElement.name,dgElement.degree,mesh);
            [x,g4e] = FormTest.elementwiseConstant(fe);
            [weight,g1,g2] = FormTest.sideIntegrals(mesh,g4e);
            inner = mesh.e4ed(:,2) > 0;

            b = LinearForm('f(x)*{v} dsi').assemble(fe,'f',@(x)1 + x(:,1));

            testCase.verifyEqual(x'*b,sum(weight(inner).*(g1(inner) + g2(inner))/2), ...
                'AbsTol',testCase.tol,'integral of f times the average');
        end

        function testParserKeywords(testCase,parserCase)
            form = BilinearForm(parserCase.str,parserCase.varStr);

            testCase.verifyEqual([form.forms.operator1,form.forms.operator2],parserCase.ops);
        end

        function testParserCoefficient(testCase)
            form = BilinearForm('mu(x)*u*v dx');

            testCase.verifyEqual(form.forms.coefficients,{'mu'});
        end

        function testCurlIdentities(testCase)
            mesh = PatchTestCase.patchMesh(struct('name','Lshape'));
            P2 = Pk(2,mesh);
            P3 = Pk(3,mesh);
            nd = Nd1k(2,mesh);

            A = BilinearForm('curl(u)*curl(v) dx').assemble(P2);
            G = BilinearForm('grad(u)*grad(v) dx').assemble(P2);
            testCase.verifyEqual(full(A),full(G),'AbsTol',testCase.tol,'|curl u| = |grad u|');
            A = BilinearForm('curlcurl(u):curlcurl(v) dx').assemble(P3);
            H = BilinearForm('d2(u):d2(v) dx').assemble(P3);
            testCase.verifyEqual(full(A),full(H),'AbsTol',testCase.tol*norm(full(H),inf), ...
                'Airy matrix and Hessian');
            free = P3.freeNodes;
            A = BilinearForm('curl(u)*v dx').assemble(P3,nd);
            B = BilinearForm('u*rot(v) dx').assemble(P3,nd);
            testCase.verifyEqual(full(A(:,free)),full(B(:,free)),'AbsTol',testCase.tol, ...
                '(curl u,v) = (u,rot v)');
        end

        function testCurlIdentity3d(testCase)
            nd = Nd1k3d(2,UnitTestCase.cubeMesh());
            free = nd.freeNodes;

            A = BilinearForm('curl(u)*v dx').assemble(nd);
            B = BilinearForm('u*curl(v) dx').assemble(nd);
            testCase.verifyEqual(full(A(free,free)),full(B(free,free)),'AbsTol',testCase.tol, ...
                '(curl u,v) = (u,curl v)');
        end

    end

    methods (Static)

        function verifySymmetric(testCase,A)
            %verifySymmetric - Verify that A is symmetric up to round-off
            testCase.verifyLessThan(norm(A-A',1),testCase.tol*norm(A,1),'symmetric');
        end

        function verifySemidefinite(testCase,A)
            %verifySemidefinite - Verify that A is positive semidefinite
            lambda = eig(full(A+A')/2);
            testCase.verifyGreaterThan(min(lambda),-testCase.tol*max(abs(lambda)), ...
                'positive semidefinite');
        end

        function verifyDefinite(testCase,A)
            %verifyDefinite - Verify that A is positive definite
            [~,flag] = chol(sparse(A+A')/2);
            testCase.verifyEqual(flag,0,'positive definite');
        end

        function [x,g4e] = elementwiseConstant(fe)
            %elementwiseConstant - Coefficients of g = 1 + T/nrElems on T
            nrElems = fe.mesh.nrElems;
            g4e = 1 + (1:nrElems)'/nrElems;
            x = fe.project(@(pts,ptsref,ind)g4e(ind(:)));
        end

        function [weight,g1,g2] = sideIntegrals(mesh,g4e)
            %sideIntegrals - Integrals of f = 1 + x_1 and side values of g
            %   WEIGHT are the integrals of f over the edges, G1 and G2 the
            %   values of g on the first and second neighboring element
            %   (zero on the boundary).
            mid = (mesh.c4n(mesh.n4ed(:,1),:) + mesh.c4n(mesh.n4ed(:,2),:))/2;
            weight = mesh.length4ed.*(1 + mid(:,1));
            e4ed = mesh.e4ed;
            g1 = g4e(e4ed(:,1));
            g2 = zeros(size(g1));
            g2(e4ed(:,2) > 0) = g4e(e4ed(e4ed(:,2) > 0,2));
        end

    end

end
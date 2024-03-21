classdef fem4mModuleEigenvalue < fem4mModuleTemplate
    %fem4mModuleEigenvalue - Adaptive FEM module for eigenvalue problems
    %   OBJ = fem4mModuleEigenvalue() creates the Eigenvalue module of the
    %   fem4m app for the Dirichlet eigenvalue problem of the Laplacian
    %   on the domains Square, Lshape and Slit, with the Lagrange
    %   elements P1 to P3, the refinement strategies 'uniform' and
    %   'adaptive', the solver 'eigs' and the tunable parameter m.
    %
    %   eigenvalue_exact holds the eigenvalues m = 1,...,10 of the Square
    %   and the first eigenvalue (m = 1) of Lshape and Slit; all other
    %   entries are 0, so evError4lvl is the eigenvalue error only for
    %   these entries.
    %
    %   fem4mModuleEigenvalue functions:
    %       loadGeometry    - Create the mesh for a selected problem
    %       afem            - Run the adaptive finite element loop
    %       drawConvergence - Plot the convergence history curves
    %       plot            - Plot the computed finite element solution
    %
    %   fem4mModuleEigenvalue properties:
    %       eigenvalue_exact - Exact eigenvalues, one row per problem
    %       m                - Index of the computed eigenpair
    %       theta            - Bulk marking parameter for adaptive refinement
    %       maxM             - Number of columns of eigenvalue_exact
    %       eigenvalue       - Computed eigenvalue
    %       eta4lvl          - Error estimator values per level
    %       nrDoF4lvl        - Number of degrees of freedom per level
    %       evError4lvl      - Eigenvalue errors per level
    %       x                - Coefficients of the computed eigenfunction
    %       fe               - Current finite element space object
    %
    %   See also fem4mModuleTemplate, fem4mModuleElliptic,
    %   fem4mModuleElasticity

    % Copyright (c) 2026 by the fem4m authors

    properties
        eigenvalue_exact % Exact eigenvalues, one row per problem
        m = 1            % Index of the computed eigenpair
        theta = 0.5;     % Bulk marking parameter for adaptive refinement
        maxM = 10;       % Number of columns of eigenvalue_exact
        eigenvalue       % Computed eigenvalue
        eta4lvl          % Error estimator values per level
        nrDoF4lvl        % Number of degrees of freedom per level
        evError4lvl      % Eigenvalue errors per level
        x                % Coefficients of the computed eigenfunction
        fe               % Current finite element space object
    end

    methods

        function obj = fem4mModuleEigenvalue()
            %fem4mModuleEigenvalue - Create the Eigenvalue app module
            obj.title = 'Eigenvalue';
            obj.problems = {
                'Square'
                'Lshape'
                'Slit'
                };
            obj.eigenvalue_exact = zeros(numel(obj.problems),obj.maxM);
            obj.eigenvalue_exact(1,:) = [2*pi^2 5*pi^2 5*pi^2 8*pi^2 10*pi^2 10*pi^2 13*pi^2 13*pi^2 17*pi^2 17*pi^2];
            obj.eigenvalue_exact(2,1) = 9.6397238440219;
            obj.eigenvalue_exact(3,1) = 8.3713297112;
            obj.finiteElements = {'P1' 'P2' 'P3'};
            obj.refine = {'uniform','adaptive'};
            obj.solver = {'eigs'};
            obj.parameter = {'m'};
        end

        function mesh = loadGeometry(obj,problemIndx)
            %loadGeometry - Create the mesh for a selected problem
            %   MESH = loadGeometry(OBJ,problemIndx) creates the mesh of
            %   the domain OBJ.problems{problemIndx} with the whole
            %   boundary (its Dirichlet and Neumann edges) as Dirichlet
            %   boundary.
            mesh = Mesh(struct('name',obj.problems{problemIndx}));
            mesh = Mesh(mesh.c4n,mesh.n4e,[mesh.n4Db;mesh.n4Nb],[]);
        end

        function afem(obj,app)
            %AFEM - Run the adaptive finite element loop
            %   AFEM(OBJ,APP) computes the OBJ.m-th smallest eigenvalue
            %   and an eigenfunction of the Laplacian with the Dirichlet
            %   boundary of APP.mesh (the whole boundary for the mesh of
            %   loadGeometry), starting on the mesh APP.mesh, with the
            %   Lagrange elements of degree APP.feIndx and eigs (shift
            %   1/5) on a sequence of meshes refined uniformly or
            %   adaptively (MarkBulk with the bulk parameter OBJ.theta and
            %   the residual error estimator) by the strategy
            %   APP.refineIndx. OBJ.m must be a positive integer of at
            %   most OBJ.maxM (10); a larger m makes afem fail when it
            %   computes evError4lvl. A mesh with fewer than OBJ.m free
            %   degrees of freedom is first refined uniformly. The loop
            %   stops after the level on which the number of degrees of
            %   freedom reaches APP.maxNrDoF or the level APP.maxLvl. For
            %   every level it stores the number of degrees of freedom,
            %   the estimator and the eigenvalue error, the distance to
            %   eigenvalue_exact(APP.problemIndx,OBJ.m), in nrDoF4lvl,
            %   eta4lvl and evError4lvl, and it keeps the eigenvalue and
            %   eigenfunction of the last level in OBJ.eigenvalue, OBJ.x
            %   and OBJ.fe.
            mesh = app.mesh;
            maxNrDoF = app.maxNrDoF;
            maxLvl = app.maxLvl;
            feIndx = app.feIndx;
            refineIndx = app.refineIndx;
            problemIndx = app.problemIndx;

            obj.nrDoF4lvl = [];
            obj.eta4lvl = [];
            obj.evError4lvl = [];
            if strcmpi(obj.refine{refineIndx},'adaptive')
                mark = MarkBulk;
                mark.theta = obj.theta;
            else
                mark = MarkUniform;
            end
            lvl = 1;
            while(true)
                while(true)
                    obj.fe = eval(['Pk(',num2str(feIndx),',mesh)']);
                    if numel(obj.fe.freeNodes)<obj.m
                        mesh = mesh.refine('uniform');
                    else
                        break
                    end
                end
                a = BilinearForm('grad(u)*grad(v)dx');
                A = a.assemble(obj.fe);
                b = BilinearForm('u*v dx');
                B = b.assemble(obj.fe);
                obj.x = zeros(size(A,1),1);
                freeNodes = obj.fe.freeNodes;
                opts = [];
                opts.disp = 0;
                [V,D] = eigs(A(freeNodes,freeNodes),B(freeNodes,freeNodes),obj.m,1/5,opts);
                [D,ind] = sort(diag(D));
                obj.x(freeNodes) = V(:,ind(obj.m));
                obj.eigenvalue = D(obj.m);
                %estimate
                eta = EstimatorForm('h^2*(f(x)+div(grad(u)))^2dx + h^1*([grad(u)*n])^2dsi');
                eta4e = eta.assemble(obj.x,obj.fe,'f',...
                    @(pts,pts_ref,inices)(obj.eigenvalue.*obj.fe.getValue([],pts_ref,[],obj.x)));
                % error
                obj.eta4lvl(end+1) = norm(eta4e);
                obj.nrDoF4lvl(end+1) = obj.fe.nrDoF;
                obj.evError4lvl(end+1) = abs(obj.eigenvalue-obj.eigenvalue_exact(problemIndx,obj.m));
                % display
                disp(['nodes/dofs: ',num2str(size(mesh.c4n,1)),'/',num2str(obj.nrDoF4lvl(end)),...
                    '; ev = ',num2str(obj.eigenvalue),...
                    '; evError = ',num2str(obj.evError4lvl(end)),...
                    '; estimator = ',num2str(obj.eta4lvl(end))]);
                if obj.plotMesh
                    app.plotMesh;
                end
                if obj.fe.nrDoF >= maxNrDoF || lvl >= maxLvl
                    break;
                end
                % mark & refine
                mark = mark.mark(obj.fe,eta4e);
                mesh = mesh.refine(mark.marked4ed);
                lvl = lvl + 1;
                app.mesh = mesh;
            end
        end

        function drawConvergence(obj,checkBox1,checkBox2,~)
            %drawConvergence - Plot the convergence history curves
            %   drawConvergence(OBJ,checkBox1,checkBox2) plots the
            %   squared estimator (checkBox1) and the eigenvalue error
            %   (checkBox2) of the last afem run against the number of
            %   degrees of freedom by plotConvergence, every curve whose
            %   logical flag is true.
            %
            %   drawConvergence(OBJ,checkBox1,checkBox2,checkBox3) does
            %   the same; the app passes a third flag, which is ignored.
            if checkBox1
                plotConvergence(obj.nrDoF4lvl,obj.eta4lvl.^2,'\eta^2');hold on;
            end
            if checkBox2
                plotConvergence(obj.nrDoF4lvl,obj.evError4lvl,'|\lambda-\lambda_h|');hold on;
            end
        end

        function plot(obj,~)
            %PLOT - Plot the computed finite element solution
            %   PLOT(OBJ) plots the eigenfunction OBJ.x of the last level
            %   on the mesh of OBJ.fe with interpolated shading. On a mesh
            %   with fewer than 500 nodes the reference element is refined
            %   min(4,K-1) times for the plot, K the degree of OBJ.fe.
            %
            %   PLOT(OBJ,problemIndx) does the same; the app passes the
            %   index problemIndx of the current problem, which is
            %   ignored.
            if obj.fe.mesh.nrNodes < 5e2
                obj.fe.plot(obj.x,min(4,obj.fe.degree-1));
            else
                obj.fe.plot(obj.x);
            end
            shading interp;
        end

    end
end
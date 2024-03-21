classdef fem4mModuleElasticity < fem4mModuleTemplate
    %fem4mModuleElasticity - Adaptive FEM module for elasticity problems
    %   OBJ = fem4mModuleElasticity() creates the Elasticity module of
    %   the fem4m app with the problems Square, Lshape, Cooks and
    %   CircularInclusion (the classes in functionName), the
    %   vector-valued Lagrange elements P1 to P4, the refinement
    %   strategies 'uniform' and 'adaptive', and the tunable material
    %   parameters E = 100000 and nu = 0.3.
    %
    %   fem4mModuleElasticity functions:
    %       loadGeometry    - Create the mesh for a selected problem
    %       afem            - Run the adaptive finite element loop
    %       drawConvergence - Plot the convergence history curves
    %       plot            - Plot the computed finite element solution
    %
    %   fem4mModuleElasticity properties:
    %       E              - Young's modulus of the material
    %       nu             - Poisson's ratio of the material
    %       x              - Coefficients of the computed solution
    %       fe             - Current vector-valued finite element space object
    %       nrDoF4lvl      - Number of degrees of freedom per level
    %       eta4lvl        - Error estimator values per level
    %       SigmaError4lvl - Stress error per level
    %       L2Error4lvl    - L2 error of the displacement per level
    %       functionName   - Cell array of problem class names
    %
    %   See also fem4mModuleTemplate, fem4mModuleElliptic,
    %   fem4mModuleEigenvalue

    % Copyright (c) 2026 by the fem4m authors

    properties
        E = 100000     % Young's modulus of the material
        nu = 0.3       % Poisson's ratio of the material
        x              % Coefficients of the computed solution
        fe             % Current vector-valued finite element space object
        nrDoF4lvl      % Number of degrees of freedom per level
        eta4lvl        % Error estimator values per level
        SigmaError4lvl % Stress error per level
        L2Error4lvl    % L2 error of the displacement per level
        functionName   % Cell array of problem class names
    end

    methods

        function obj = fem4mModuleElasticity()
            %fem4mModuleElasticity - Create the Elasticity app module
            obj.title = 'Elasticity';
            obj.problems = {
                'Square'
                'Lshape'
                'Cooks'
                'CircularInclusion'
                };
            obj.functionName = {
                'problemElasticitySquareExactLocking'
                'problemElasticityLshapeExact'
                'problemElasticityCooks'
                'problemElasticityCircularInclusion'
                };
            obj.finiteElements = {'P1' 'P2' 'P3' 'P4'};
            obj.refine = {'uniform','adaptive'};
            obj.parameter = {'E','nu'};
        end

        function mesh = loadGeometry(obj,problemIndx)
            %loadGeometry - Create the mesh for a selected problem
            %   MESH = loadGeometry(OBJ,problemIndx) creates an object of
            %   the problem class OBJ.functionName{problemIndx} with the
            %   material parameters OBJ.E and OBJ.nu and returns the mesh
            %   of its geometry.
            p = eval([obj.functionName{problemIndx},'(''E'',obj.E,''nu'',obj.nu)']);
            mesh = Mesh(p.geom);
        end

        function afem(obj,app)
            %AFEM - Run the adaptive finite element loop
            %   AFEM(OBJ,APP) solves the linear elasticity problem
            %   selected by APP.problemIndx with the material parameters
            %   OBJ.E and OBJ.nu, starting on the mesh APP.mesh, with the
            %   vector-valued Lagrange elements (PkVector) of degree
            %   APP.feIndx on a sequence of meshes refined uniformly or
            %   adaptively (MarkBulk with the residual error estimator) by
            %   the strategy APP.refineIndx. The loop stops after the
            %   level on which the number of degrees of freedom reaches
            %   APP.maxNrDoF or the level APP.maxLvl. For every level it
            %   stores the number of degrees of freedom and the estimator
            %   in nrDoF4lvl and eta4lvl and, if the problem has an exact
            %   solution (a method exactU), the L2 error of the
            %   displacement and the error of the stress C*eps(u_h) in the
            %   L2 norm weighted by the inverse of C in L2Error4lvl and
            %   SigmaError4lvl; otherwise (Cooks) these stay empty. It
            %   keeps the solution of the last level in OBJ.x and OBJ.fe.
            mesh = app.mesh;
            maxNrDoF = app.maxNrDoF;
            maxLvl = app.maxLvl;
            feIndx = app.feIndx;
            refineIndx = app.refineIndx;
            problemIndx = app.problemIndx;

            obj.nrDoF4lvl = [];
            obj.eta4lvl = [];
            obj.L2Error4lvl = [];
            obj.SigmaError4lvl = [];
            if strcmpi(obj.refine{refineIndx},'adaptive')
                mark = MarkBulk;
            else
                mark = MarkUniform;
            end
            p = eval([obj.functionName{problemIndx},'(''E'',obj.E,''nu'',obj.nu)']);
            lvl = 1;
            while(true)
                obj.fe = eval(['PkVector(',num2str(feIndx),',mesh)']);
                % assembly
                a = BilinearForm('C(x)*eps(u):eps(v)dx');
                A = a.assemble(obj.fe,'C',p.C);
                b = LinearForm('f(x)*v dx + g(x)*vdsn');
                normals = mesh.normal4ed;
                b = b.assemble(obj.fe,'f',@p.f,'g',@(x)(p.g(x,normals(mesh.ed4Nb,:))));

                u4Db = obj.fe.interpolate(@p.u4Db);
                obj.x = zeros(obj.fe.nrDoF,1);
                obj.x(obj.fe.fixedNodes) = u4Db(obj.fe.fixedNodes);
                b = b - A*obj.x;
                % solve
                freeNodes = obj.fe.freeNodes;
                obj.x(freeNodes) = A(freeNodes,freeNodes) \ b(freeNodes);
                % estimate
                eta = EstimatorForm(['h^2*(f(x)+div(C(x)eps(u))^2dx', ...
                    '+ h^1*([C(x)*eps(u)*n])^2dsi', ...
                    '+ h^1*(g(x)-C(x)*eps(u)*n)^2dsn']);
                eta4e = eta.assemble(obj.x,obj.fe,'f',@p.f,'g',@(x,n)p.g(x,n),'C',p.C);
                % error
                obj.eta4lvl(end+1) = norm(eta4e);
                obj.nrDoF4lvl(end+1) = obj.fe.nrDoF;
                if ismethod(p,'exactU')
                    error = ErrorForm('(exact(x)-u)^2 dx');
                    obj.L2Error4lvl(end+1) = error.assemble(obj.x,obj.fe,2*obj.fe.degree+2,'exact',@p.exactU);
                    error = ErrorForm('(exact(x)-C(x)*eps(u))^2Cin(x)dx');
                    obj.SigmaError4lvl(end+1) =  error.assemble(obj.x,obj.fe,2*obj.fe.degree+2,'C',p.C,'exact',@p.exactSigma,'Cin',p.CInv);
                    disp(['nodes/dofs: ',num2str(size(mesh.c4n,1)),'/',num2str(obj.nrDoF4lvl(end)),...
                        '; DisplacementError = ',num2str(obj.L2Error4lvl(end)),...
                        '; SigmaError = ',num2str(obj.SigmaError4lvl(end)),...
                        '; Estimator = ',num2str(obj.eta4lvl(end))]);
                else
                    disp(['nodes/dofs: ',num2str(size(mesh.c4n,1)),'/',num2str(obj.nrDoF4lvl(end)),...
                        '; Estimator = ',num2str(obj.eta4lvl(end))]);
                end
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

        function drawConvergence(obj,checkBox1,checkBox2,checkBox3)
            %drawConvergence - Plot the convergence history curves
            %   drawConvergence(OBJ,checkBox1,checkBox2,checkBox3) plots
            %   the estimator (checkBox1), the stress error (checkBox2)
            %   and the displacement error (checkBox3) of the last afem
            %   run against the number of degrees of freedom by
            %   plotConvergence, every curve whose logical flag is true;
            %   the error curves only if the errors were computed.
            if checkBox1
                plotConvergence(obj.nrDoF4lvl,obj.eta4lvl,'\eta');hold on;
            end
            if checkBox2 && ~isempty(obj.SigmaError4lvl)
                plotConvergence(obj.nrDoF4lvl,obj.SigmaError4lvl,'||\nabla(u-u_h)||');hold on;
            end
            if checkBox3 && ~isempty(obj.L2Error4lvl)
                plotConvergence(obj.nrDoF4lvl,obj.L2Error4lvl,'||u-u_h||');hold on;
            end
        end

        function plot(obj,problemIndx)
            %PLOT - Plot the computed finite element solution
            %   PLOT(OBJ,problemIndx) plots the mesh of OBJ.fe displaced by
            %   the last computed displacement OBJ.x, amplified by the
            %   factor scale of the problem selected by problemIndx, and
            %   colored by the von Mises stress of C*eps(u_h) interpolated
            %   into DGPkSymMatrix of the degree of OBJ.fe, with
            %   interpolated shading. problemIndx must select the problem
            %   of the last afem run. For every problem but the first, the
            %   color limits are clipped to [0,5].
            p = eval([obj.functionName{problemIndx},'(''E'',obj.E,''nu'',obj.nu)']);
            pkSigma = DGPkSymMatrix(obj.fe.degree,obj.fe.mesh);
            sigma = pkSigma.interpolate(@(pts,pts_ref,parts)matMul(permute(p.C,[3 1 2]),obj.fe.getSymGradValue(pts,pts_ref,parts,obj.x)));
            obj.fe.plotDisplacement(obj.x,p.scale,sigma,pkSigma);
            if problemIndx ~= 1
                limit = clim;
                clim([max(0,limit(1)),min(5,limit(2))]);
            end
            shading interp;
        end

    end
end
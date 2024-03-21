classdef fem4mModuleElliptic < fem4mModuleTemplate
    %fem4mModuleElliptic - Adaptive FEM module for scalar elliptic problems
    %   OBJ = fem4mModuleElliptic() creates the Elliptic module of the
    %   fem4m app with the problems Square, Lshape, Slit,
    %   SquareWaterfall, SquareFullElliptic, Kellogg and MountainValley
    %   (the classes in functionName), the Lagrange elements P1 to P4,
    %   the refinement strategies 'uniform', 'adaptive' and 'graded',
    %   and no tunable parameters.
    %
    %   fem4mModuleElliptic functions:
    %       loadGeometry    - Create the mesh for a selected problem
    %       afem            - Run the adaptive finite element loop
    %       drawConvergence - Plot the convergence history curves
    %       plot            - Plot the computed finite element solution
    %
    %   fem4mModuleElliptic properties:
    %       x              - Coefficients of the computed solution
    %       fe             - Current finite element space object
    %       nrDoF4lvl      - Number of degrees of freedom per level
    %       eta4lvl        - Error estimator values per level
    %       SigmaError4lvl - L2 error of the gradient per level
    %       L2Error4lvl    - L2 error of the solution per level
    %       functionName   - Cell array of problem class names
    %
    %   See also fem4mModuleTemplate, fem4mModuleElasticity,
    %   fem4mModuleEigenvalue

    % Copyright (c) 2026 by the fem4m authors

    properties
        x              % Coefficients of the computed solution
        fe             % Current finite element space object
        nrDoF4lvl      % Number of degrees of freedom per level
        eta4lvl        % Error estimator values per level
        SigmaError4lvl % L2 error of the gradient per level
        L2Error4lvl    % L2 error of the solution per level
        functionName   % Cell array of problem class names
    end

    methods

        function obj = fem4mModuleElliptic()
            %fem4mModuleElliptic - Create the Elliptic app module
            obj.title = 'Elliptic';
            obj.problems = {
                'Square'
                'Lshape'
                'Slit'
                'SquareWaterfall'
                'SquareFullElliptic'
                'Kellogg'
                'MountainValley'
                };
            obj.functionName = {
                'problemScalarEllipticSquareExact'
                'problemScalarEllipticLshapeExact'
                'problemScalarEllipticSlitExact'
                'problemScalarEllipticSquareWaterfallExact'
                'problemScalarEllipticSquareFullEllipticExact'
                'problemScalarEllipticSquareKelloggExact'
                'problemScalarEllipticSquareMountainValleyExact'
                };
            obj.finiteElements = {'P1' 'P2' 'P3' 'P4'};
            obj.refine = {'uniform','adaptive','graded'};
            obj.parameter = {};
        end

        function mesh = loadGeometry(obj,problemIndx)
            %loadGeometry - Create the mesh for a selected problem
            %   MESH = loadGeometry(OBJ,problemIndx) creates an object of
            %   the problem class OBJ.functionName{problemIndx} and
            %   returns the mesh of its geometry, refined once uniformly
            %   for the SquareFullElliptic and Kellogg problems.
            p = eval(obj.functionName{problemIndx});
            mesh = Mesh(p.geom);
            if strcmpi(obj.functionName{problemIndx},'problemScalarEllipticSquareFullEllipticExact') || ...
                    strcmpi(obj.functionName{problemIndx},'problemScalarEllipticSquareKelloggExact')
                mesh = mesh.refine('uniform');
            end
        end

        function afem(obj,app)
            %AFEM - Run the adaptive finite element loop
            %   AFEM(OBJ,APP) solves the scalar elliptic problem selected
            %   by APP.problemIndx, starting on the mesh APP.mesh, with
            %   the Lagrange elements of degree APP.feIndx on a sequence
            %   of meshes refined by the strategy APP.refineIndx:
            %   'uniform' refines every element, 'adaptive' marks by
            %   MarkBulk with the residual error estimator, and 'graded'
            %   refines uniformly and, for Lshape and Slit, moves the
            %   nodes within the radius 1/2 of the origin towards it
            %   (refine with 'gradedbeta', mu = 2/3 for Lshape and
            %   mu = 1/4 for Slit); for every other problem 'graded'
            %   refines uniformly. The loop stops after the level on which
            %   the number of degrees of freedom reaches APP.maxNrDoF or
            %   the level APP.maxLvl. For every level it stores the
            %   number of degrees of freedom, the estimator, the L2 error
            %   and the L2 error of the gradient (H1 seminorm of the error)
            %   in nrDoF4lvl, eta4lvl, L2Error4lvl and SigmaError4lvl, and
            %   it keeps the solution of the last level in OBJ.x and OBJ.fe.
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
            p = eval(obj.functionName{problemIndx});
            lvl = 1;
            while(true)
                obj.fe = eval(['Pk(',num2str(feIndx),',mesh)']);
                % solve
                a = BilinearForm('alpha(x)grad(u)*grad(v)dx + beta(x)*grad(u)vdx + gamma(x)u*v dx');
                A = a.assemble(obj.fe,'alpha',@p.alpha,'beta',@p.beta,'gamma',@p.gamma);
                rhs = LinearForm('f(x)vdx + g(x)vdsn');
                b = rhs.assemble(obj.fe,'f',@p.f,'g',@(x,n)p.g(x,n));
                obj.x = zeros(size(A,1),1);
                c4dof = obj.fe.c4dof;
                obj.x(obj.fe.fixedNodes) = p.u4Db(c4dof(obj.fe.fixedNodes,:));
                b = b - A*obj.x;
                obj.x(obj.fe.freeNodes) = A(obj.fe.freeNodes,obj.fe.freeNodes)\b(obj.fe.freeNodes);
                % estimate
                eta = EstimatorForm(['h^2*p^-2*(f(x)+div(alphac(x)grad(u))-beta(x)*grad(u)-gamma(x)u)^2dx', ...
                    '+ h^1*p^-1*([alpha(x)*grad(u)*n])^2dsi', ...
                    '+ h^1*p^-1*(g(x)-alpha(x)*grad(u)*n)^2dsn']);
                eta4e = eta.assemble(obj.x,obj.fe,'f',@p.f,'g',@(x,n)p.g(x,n),...
                    'alphac',{@p.alpha,@p.dalpha},'alpha',@p.alpha,'beta',@p.beta,'gamma',@p.gamma);
                % error
                obj.eta4lvl(end+1) = norm(eta4e);
                obj.nrDoF4lvl(end+1) = obj.fe.nrDoF;
                l2 = ErrorForm('(exact(x)-u)^2dx');
                obj.L2Error4lvl(end+1) = l2.assemble(obj.x,obj.fe,2*obj.fe.degree+2,'exact',@p.exactU);
                h1 = ErrorForm('(exactGrad(x)-grad(u))^2dx');
                obj.SigmaError4lvl(end+1) = h1.assemble(obj.x,obj.fe,2*obj.fe.degree,'exactGrad',@p.exactGradU);
                % display
                disp(['nodes/dofs: ',num2str(size(mesh.c4n,1)),'/',num2str(obj.nrDoF4lvl(end)),...
                    '; L2Error = ',num2str(obj.L2Error4lvl(end)),...
                    '; H1Error = ',num2str(obj.SigmaError4lvl(end)),...
                    '; estimator = ',num2str(obj.eta4lvl(end))]);
                if obj.plotMesh
                    app.plotMesh;
                end
                if obj.fe.nrDoF >= maxNrDoF || lvl >= maxLvl
                    break;
                end
                % mark & refine
                if strcmpi(obj.refine{refineIndx},'graded')
                    if strcmpi(obj.functionName{problemIndx},'problemScalarEllipticLshapeExact')
                        mesh = mesh.refine('gradedbeta',2/3,[0,0],1/2);
                    elseif strcmpi(obj.functionName{problemIndx},'problemScalarEllipticSlitExact')
                        mesh = mesh.refine('gradedbeta',1/4,[0,0],1/2);
                    else
                        mesh = mesh.refine('uniform');
                    end
                else
                    mark = mark.mark(obj.fe,eta4e);
                    mesh = mesh.refine(mark.marked4ed);
                end
                lvl = lvl + 1;
                app.mesh = mesh;
            end
        end

        function drawConvergence(obj,checkBox1,checkBox2,checkBox3)
            %drawConvergence - Plot the convergence history curves
            %   drawConvergence(OBJ,checkBox1,checkBox2,checkBox3) plots
            %   the estimator (checkBox1), the L2 error of the gradient
            %   (checkBox2) and the L2 error (checkBox3) of the last afem
            %   run against the number of degrees of freedom by
            %   plotConvergence, every curve whose logical flag is true.
            if checkBox1
                plotConvergence(obj.nrDoF4lvl,obj.eta4lvl,'\eta');hold on;
            end
            if checkBox2
                plotConvergence(obj.nrDoF4lvl,obj.SigmaError4lvl,'||\nabla(u-u_h)||');hold on;
            end
            if checkBox3
                plotConvergence(obj.nrDoF4lvl,obj.L2Error4lvl,'||u-u_h||');hold on;
            end
        end

        function plot(obj,~)
            %PLOT - Plot the computed finite element solution
            %   PLOT(OBJ) plots the solution OBJ.x of the last level on
            %   the mesh of OBJ.fe with interpolated shading. On a mesh
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
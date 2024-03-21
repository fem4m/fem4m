classdef fem4mModuleTemplate < handle
    %fem4mModuleTemplate - Base class for fem4m app problem modules
    %   fem4mModuleTemplate is an abstract handle class. A module of the
    %   fem4m app subclasses it, sets the properties in its constructor
    %   and implements the abstract methods loadGeometry, afem,
    %   drawConvergence and plot. The app fem4m.mlapp creates one object
    %   of every class fem4mModule*.m in this folder except
    %   fem4mModuleTemplate, by its constructor without arguments.
    %
    %   The argument APP of afem and save is the app or a struct with the
    %   fields mesh (the current mesh), maxNrDoF (limit of the degrees of
    %   freedom), maxLvl (limit of the levels), and feIndx, refineIndx,
    %   problemIndx, estimateIndx and solverIndx (indices into
    %   finiteElements, refine, problems, estimate and solver). When
    %   OBJ.plotMesh is true, afem also calls APP.plotMesh, a method that
    %   plots APP.mesh, after every level. afem assigns every refined
    %   mesh to APP.mesh, which the caller sees only if APP is a handle
    %   object such as the app. The app calls drawConvergence with three
    %   logical flags that select the convergence curves, the estimator
    %   first, and plot with the index problemIndx of the current
    %   problem as second argument.
    %
    %   save stores every nonempty numeric property whose name ends in
    %   4lvl, so a module keeps its level-wise results in such properties
    %   and must have the property nrDoF4lvl.
    %
    %   fem4mModuleTemplate functions:
    %       loadGeometry    - Create the mesh for a selected problem
    %       afem            - Run the adaptive finite element loop
    %       drawConvergence - Plot the convergence history curves
    %       plot            - Plot the computed finite element solution
    %       result          - Results file description of the settings
    %       getFileName     - Build the results file name from settings
    %       save            - Save the level-wise results to a MAT-file
    %       batchRun        - Run all experiments in serial or parallel
    %       runAndSave      - Run one experiment and save its results
    %
    %   fem4mModuleTemplate properties:
    %       title          - Title string of the module
    %       problems       - Cell array of problem names
    %       finiteElements - Cell array of finite element space names
    %       refine         - Cell array of refinement strategies
    %       estimate       - Cell array of error estimator names
    %       solver         - Cell array of solver names
    %       parameter      - Cell array of tunable parameter property names
    %       plotMesh       - Flag to plot the mesh during refinement
    %
    %   See also fem4mModuleElliptic, fem4mModuleElasticity,
    %   fem4mModuleEigenvalue, fem4mResult

    % Copyright (c) 2026 by the fem4m authors

    properties
        title          % Title string of the module
        problems       % Cell array of problem names
        finiteElements % Cell array of finite element space names
        refine = {'uniform'}     % Cell array of refinement strategies
        estimate = {'residual'}  % Cell array of error estimator names
        solver = {'direct'}      % Cell array of solver names
        parameter      % Cell array of tunable parameter property names
        plotMesh = true          % Flag to plot the mesh during refinement
    end

    methods (Abstract)

        %loadGeometry - Create the mesh for a selected problem
        mesh = loadGeometry(obj,problemIndx);
        %AFEM - Run the adaptive finite element loop
        afem(obj,app);
        %drawConvergence - Plot the convergence history curves
        drawConvergence(obj,checkBox1,checkBox2,checkBox3);
        %PLOT - Plot the computed finite element solution
        plot(obj);

    end

    methods

        function r = result(obj,maxNrDoF,feIndx,refineIndx,problemIndx,estimateIndx,solverIndx)
            %RESULT - Results file description of the settings
            %   R = RESULT(OBJ,maxNrDoF,feIndx,refineIndx,problemIndx,estimateIndx,solverIndx)
            %   returns the fem4mResult of the problem, finite element
            %   space, refinement strategy, error estimator and solver
            %   selected by the indices into OBJ.problems,
            %   OBJ.finiteElements, OBJ.refine, OBJ.estimate and
            %   OBJ.solver, with the limit maxNrDoF. The values of the
            %   properties named by the first two entries of
            %   OBJ.parameter give param1 and param2, their names the
            %   labels; a missing entry leaves the slot 0 with an empty
            %   label.
            r = fem4mResult(problem=obj.problems{problemIndx}, ...
                fe=obj.finiteElements{feIndx}, ...
                refine=obj.refine{refineIndx}, ...
                estimator=obj.estimate{estimateIndx}, ...
                solver=obj.solver{solverIndx}, ...
                param1=0,param2=0,maxNrDoF=maxNrDoF);
            if numel(obj.parameter) > 0
                r.param1 = obj.(obj.parameter{1});
                r.label.param1 = obj.parameter{1};
            end
            if numel(obj.parameter) > 1
                r.param2 = obj.(obj.parameter{2});
                r.label.param2 = obj.parameter{2};
            end
        end

        function val = getFileName(obj,maxNrDoF,feIndx,refineIndx,problemIndx,estimateIndx,solverIndx)
            %getFileName - Build the results file name from settings
            %   VAL = getFileName(OBJ,maxNrDoF,feIndx,refineIndx,problemIndx,estimateIndx,solverIndx)
            %   returns the file name, without extension, of the results
            %   of the selected settings: fileName of the fem4mResult
            %   that result returns for the same arguments.
            val = obj.result(maxNrDoF,feIndx,refineIndx,problemIndx, ...
                estimateIndx,solverIndx).fileName();
        end

        function save(obj,app)
            %SAVE - Save the level-wise results to a MAT-file
            %   SAVE(OBJ,APP) collects the numeric properties of OBJ whose
            %   names end in 4lvl into a structure, with the suffix 4lvl
            %   removed from the field names, and writes it by
            %   fem4mResult.save to the subfolder OBJ.title of the results
            %   folder (fem4mPath.resultFolder). The file name is that of
            %   getFileName for the fields maxNrDoF, feIndx, refineIndx,
            %   problemIndx, estimateIndx and solverIndx of APP. OBJ must
            %   have a nonempty nrDoF4lvl, and all saved arrays must have
            %   one entry per level. Empty arrays, e.g. the errors of a
            %   problem without exact solution, are not saved.
            prop = properties(obj);
            levels = struct();
            for j = 1 : numel(prop)
                if endsWith(prop{j},'4lvl') && isnumeric(obj.(prop{j})) ...
                        && ~isempty(obj.(prop{j}))
                    levels.(erase(prop{j},'4lvl')) = obj.(prop{j});
                end
            end

            r = obj.result(app.maxNrDoF,app.feIndx,app.refineIndx, ...
                app.problemIndx,app.estimateIndx,app.solverIndx);
            r.save(fem4mPath.resultFolder(obj.title),levels);
        end

        function batchRun(obj,maxNrDoF,ntasks,nthreads)
            %batchRun - Run all experiments in serial or parallel
            %   batchRun(OBJ) runs, one after the other, by runAndSave
            %   with maxNrDoF = 1e6, every experiment of OBJ, i.e. every
            %   combination of problem, finite element space, refinement
            %   strategy, estimator and solver with the current values of
            %   the parameters named in OBJ.parameter, whose results file
            %   (see getFileName; the names are compared ignoring case)
            %   is not yet in the subfolder OBJ.title of the results
            %   folder.
            %
            %   batchRun(OBJ,maxNrDoF) uses the limit maxNrDoF of the
            %   degrees of freedom per experiment. An empty maxNrDoF,
            %   NTASKS or NTHREADS takes its default.
            %
            %   batchRun(OBJ,maxNrDoF,NTASKS) runs the experiments by
            %   parfeval on a local parallel pool of min(NTASKS,number of
            %   remaining experiments) workers if NTASKS > 1 and displays
            %   the output of every finished experiment; NTASKS <= 1
            %   (default 1) runs them one after the other. No parallel pool may be open
            %   before; the pool is deleted at the end and on an error.
            %
            %   batchRun(OBJ,maxNrDoF,NTASKS,NTHREADS) sets the number of
            %   computational threads to NTHREADS (default 4) on every
            %   worker and on the client MATLAB session (by pctRunOnAll,
            %   not restored afterwards) when NTASKS > 1.

            if nargin < 2 || isempty(maxNrDoF)
                maxNrDoF = 1e6;
            end

            if nargin < 3 || isempty(ntasks)
                ntasks = 1;
            end

            if nargin < 4 || isempty(nthreads)
                nthreads = 4;
            end

            path = fem4mPath.resultFolder(obj.title);

            problemIndx = 1 : numel(obj.problems);
            feIndx = 1 : numel(obj.finiteElements);
            refineIndx = 1 : numel(obj.refine);
            estimateIndx = 1 : numel(obj.estimate);
            solverIndx = 1 : numel(obj.solver);

            [feIndx,refineIndx,problemIndx,estimateIndx,solverIndx] = ...
                getBatchMatrix(feIndx,refineIndx,problemIndx,estimateIndx,solverIndx);

            indx = getBatchMatrixIndicesToCompute(path,numel(feIndx),...
                @(k)obj.getFileName(maxNrDoF,feIndx(k),refineIndx(k),problemIndx(k),...
                estimateIndx(k),solverIndx(k)));
            disp(['compute remaining ', num2str(numel(indx)), ' experiments.']);
            disp('%%%%%%%%%%%%%%%% BEGIN BATCH RUN %%%%%%%%%%%%%%%%%%%');
            if ntasks <= 1
                for k = 1:numel(indx)
                    j = indx(k);
                    obj.runAndSave(maxNrDoF,feIndx(j),...
                        refineIndx(j),problemIndx(j),estimateIndx(j),solverIndx(j));
                end
            elseif numel(indx)>0
                try
                    parpool('local',max(1,min(numel(indx),ntasks)));
                    pctRunOnAll(['maxNumCompThreads(',num2str(nthreads),');']);

                    for k = 1:numel(indx)
                        j = indx(k);
                        f(k) = parfeval(@obj.runAndSave,0,maxNrDoF,feIndx(j),...
                            refineIndx(j),problemIndx(j),estimateIndx(j),solverIndx(j)); %#ok<AGROW>
                    end

                    for indx = 1 : numel(f)
                        completed = fetchNext(f);
                        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
                        disp(['computation ',num2str(completed),' finished']);
                        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
                        disp(f(completed).Diary);
                    end
                catch ME
                    delete(gcp('nocreate'));
                    error(ME.message);
                end
                delete(gcp('nocreate'));
            end
            disp('%%%%%%%%%%%%%%%% END BATCH RUN %%%%%%%%%%%%%%%%%%%');

        end

        function runAndSave(obj,maxNrDoF,feIndx,refineIndx,problemIndx,estimateIndx,solverIndx)
            %runAndSave - Run one experiment and save its results
            %   runAndSave(OBJ,maxNrDoF,feIndx,refineIndx,problemIndx,estimateIndx,solverIndx)
            %   creates the mesh of the selected problem by loadGeometry,
            %   runs afem with the selected settings and no level limit
            %   (maxLvl = inf) until the number of degrees of freedom
            %   reaches maxNrDoF, without plotting the mesh (OBJ.plotMesh
            %   is restored after a successful run), and saves the
            %   results by save.
            app = struct('maxNrDoF',maxNrDoF,'maxLvl',inf,'feIndx',feIndx, ...
                'refineIndx',refineIndx,'problemIndx',problemIndx, ...
                'estimateIndx',estimateIndx,'solverIndx',solverIndx);
            app.mesh = obj.loadGeometry(problemIndx);
            plotMesh = obj.plotMesh;
            obj.plotMesh = false;
            obj.afem(app);
            obj.plotMesh = plotMesh;
            obj.save(app);
        end

    end

end
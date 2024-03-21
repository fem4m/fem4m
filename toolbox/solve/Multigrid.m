classdef Multigrid
    %MULTIGRID - Geometric multigrid solver
    %   OBJ = MULTIGRID() creates a geometric multigrid solver without
    %   levels that uses the Richardson smoother with damping factor 1,
    %   2 pre- and 2 post-smoothing steps per level, no limit on the
    %   number of levels, and 1 coarse grid correction per level and
    %   cycle (V-cycle).
    %
    %   OBJ = MULTIGRID(SMOOTHER) uses the smoother SMOOTHER, either
    %   'Richardson', 'Jacobi', or a function handle (see
    %   setupSmoother). The damping factor then defaults to 1/8 for
    %   'Richardson', 1/2 for 'Jacobi', and 1 for a function handle or
    %   an empty SMOOTHER, which selects 'Richardson'.
    %
    %   OBJ = MULTIGRID(SMOOTHER,DAMP) also specifies the default
    %   damping factor DAMP of the levels.
    %
    %   OBJ = MULTIGRID(...,M) also specifies the default number M of
    %   pre- and post-smoothing steps per level (default 2).
    %
    %   OBJ = MULTIGRID(...,maxLvl) also keeps at most maxLvl levels
    %   (default Inf): adding a further level drops the coarsest one.
    %
    %   OBJ = MULTIGRID(...,W) also specifies the number W of coarse
    %   grid corrections per level and cycle (default 1): 1 gives a
    %   V-cycle, 2 a W-cycle. Empty arguments select the defaults.
    %
    %   Add the levels with addLevel from the coarsest to the finest
    %   mesh, then solve solves the system of the finest level. The
    %   matrices and prolongations are those of the free degrees of
    %   freedom, e.g. A(fe.freeNodes,fe.freeNodes) and
    %   P(fe.freeNodes,feOld.freeNodes) with P = fe.prolongation(feOld),
    %   or full ones together with the mask fixedNodes of addLevel.
    %
    %   MULTIGRID functions:
    %       addLevel           - Add a grid level
    %       setupSmoother      - Select the smoothing method
    %       setupLevelMatrices - Set up the fixed entries and the coarse LU
    %       contractionNumber  - Estimate the contraction number
    %       dampingFactor      - Estimate the spectral radius of T*A
    %       solve              - Solve with the multigrid iteration
    %       MG                 - Recursive multigrid cycle
    %       Richardson         - Richardson smoother
    %       Jacobi             - Jacobi smoother
    %       Poisson            - Cycle of a second multigrid solver as smoother
    %
    %   MULTIGRID properties:
    %       A          - System matrices, one per level
    %       prolong    - Prolongation operators between levels
    %       L          - Lower factor of the coarsest-level LU
    %       U          - Upper factor of the coarsest-level LU
    %       P          - Row permutation of the coarsest-level LU
    %       Q          - Column permutation of the coarsest-level LU
    %       fixedNodes - Fixed-node masks, one per level
    %       m          - Number of smoothing steps per level
    %       constM     - Default number of smoothing steps
    %       w          - Number of coarse grid corrections
    %       smoother   - Handle to the smoothing method
    %       maxLvl     - Maximum number of levels
    %       damp       - Damping factors, one per level
    %       constDamp  - Default damping factor
    %
    %   See also Schwarz, pinvit, pcg
    %
    %   Note: SMOOTHER must be a character vector or a function handle;
    %   a string scalar such as "Jacobi" raises an error.

    % Copyright (c) 2026 by the fem4m authors

    properties
        A           % System matrices, one per level
        prolong     % Prolongation operators between levels
        L           % Lower factor of the coarsest-level LU
        U           % Upper factor of the coarsest-level LU
        P           % Row permutation of the coarsest-level LU
        Q           % Column permutation of the coarsest-level LU
        fixedNodes  % Fixed-node masks, one per level
        m           % Number of smoothing steps per level
        constM      % Default number of smoothing steps
        w           % Number of coarse grid corrections
        smoother    % Handle to the smoothing method
        maxLvl      % Maximum number of levels
        damp        % Damping factors, one per level
        constDamp   % Default damping factor
    end

    methods

        function obj = Multigrid(smoother,damp,m,maxLvl,w)
            %MULTIGRID - Create multigrid solver
            if nargin < 1
                smoother = [];
            end
            if nargin < 2 || isempty(damp)
                if strcmpi(smoother,'Richardson')
                    damp = 1/8;
                elseif strcmpi(smoother,'Jacobi')
                    damp = 1/2;
                else
                    damp = 1;
                end
            end
            obj.constDamp = damp;
            if nargin < 3 || isempty(m)
                m = 2;
            end
            obj.constM = m;
            if nargin < 4 || isempty(maxLvl)
                maxLvl = inf;
            end
            obj.maxLvl = maxLvl;
            if nargin < 5 || isempty(w)
                w = 1;
            end
            obj.w = w;
            obj = obj.setupSmoother(smoother);
        end

        function obj = addLevel(obj,A,prolong,damp,m,fixedNodes)
            %addLevel - Add a grid level
            %   OBJ = addLevel(OBJ,A) adds the first (coarsest) level with the
            %   system matrix A.
            %
            %   OBJ = addLevel(OBJ,A,PROLONG) adds a finer level with the
            %   matrix A and the prolongation PROLONG from the previously
            %   added level to this one, a size(A,1)-by-N matrix with N the
            %   size of the previous matrix. PROLONG is required for every
            %   level but the first; for the first level it is ignored (pass
            %   []). When the solver already has maxLvl levels, the coarsest
            %   one is dropped.
            %
            %   OBJ = addLevel(OBJ,A,PROLONG,DAMP) also specifies the damping
            %   factor of the level, a scalar or a vector with one entry per
            %   row of A (default: the DAMP of the constructor).
            %
            %   OBJ = addLevel(OBJ,A,...,M) also specifies the number of pre-
            %   and post-smoothing steps of the level (default: the M of the
            %   constructor).
            %
            %   OBJ = addLevel(OBJ,A,...,fixedNodes) also specifies the logical
            %   mask fixedNodes (size(A,1)-by-1) of the fixed entries of the new
            %   level (default: none). The masks of the coarser levels are
            %   derived from it (a coarse entry is fixed when the sum of its
            %   prolongation weights on fixed entries is positive, i.e.
            %   (PROLONG'*fixedNodes) > 0); masks given with earlier levels
            %   are discarded. Empty arguments select the defaults.
            if isempty(A)
                error('fem4m:Multigrid:emptyMatrix','Multigrid: Add level A is empty');
            end
            computeCoarseLU = isempty(obj.L);
            if numel(obj.A) == obj.maxLvl
                computeCoarseLU = true;
                for curLvl = 1 : obj.maxLvl - 1
                    obj.A{curLvl} = obj.A{curLvl+1};
                    obj.damp{curLvl} = obj.damp{curLvl+1};
                    obj.m{curLvl} = obj.m{curLvl+1};
                end
                for curLvl = 1 : obj.maxLvl - 2
                    obj.prolong{curLvl} = obj.prolong{curLvl+1};
                end
                obj.A(end) = [];
                obj.prolong(end) = [];
                obj.damp(end) = [];
                obj.m(end) = [];
            end
            obj.A{end+1} = A;
            if numel(obj.A) > 1
                obj.prolong{end+1} = prolong;
            end
            if nargin < 4 || isempty(damp)
                damp = obj.constDamp;
            end
            if isscalar(damp)
                damp = damp*ones(size(A,1),1);
            end
            obj.damp{end+1} = damp;
            if nargin < 5 || isempty(m)
                m = obj.constM;
            end
            obj.m{end+1} = m;
            if nargin < 6 || isempty(fixedNodes)
                fixedNodes = false(size(A,1),1);
            end
            obj = obj.setupLevelMatrices(fixedNodes,computeCoarseLU);
        end

        function obj = setupSmoother(obj,smoother)
            %setupSmoother - Select the smoothing method
            %   OBJ = setupSmoother(OBJ) selects the Richardson smoother.
            %
            %   OBJ = setupSmoother(OBJ,SMOOTHER) sets the smoother to SMOOTHER,
            %   either 'Richardson' or 'Jacobi' (character vectors), or a
            %   function handle; empty selects 'Richardson'. The handle is
            %   called on each level but the coarsest (which is solved with
            %   the LU factors) as X = SMOOTHER(X,B,A,DAMP,FREE,LEVEL) with
            %   the matrix A and the damping factors DAMP of the level
            %   restricted to its free entries, the vectors X and B on these
            %   entries, the logical mask FREE of the free entries, and the
            %   level index LEVEL (1 is the coarsest), and returns the smoothed
            %   iterate. The damping factors are not changed.
            if nargin < 2 || isempty(smoother)
                obj.smoother = @Multigrid.Richardson;
            elseif isa(smoother,'function_handle')
                obj.smoother = smoother;
            elseif isa(smoother,'char')
                if strcmpi(smoother,'Richardson')
                    obj.smoother = @Multigrid.Richardson;
                elseif strcmpi(smoother,'Jacobi')
                    obj.smoother = @Multigrid.Jacobi;
                else
                    error('fem4m:Multigrid:unknownSmoother','Multigrid: Smoother (%s) not implemented',smoother);
                end
            else
                error('fem4m:Multigrid:invalidSmoother','Multigrid: parse error for variable smoother');
            end
        end

        function obj = setupLevelMatrices(obj,fixedNodes,computeCoarseLU)
            %setupLevelMatrices - Set up the fixed entries and the coarse LU
            %   OBJ = setupLevelMatrices(OBJ,fixedNodes) sets the logical mask
            %   fixedNodes of the fixed entries of the finest level, derives
            %   the masks of the coarser levels with the prolongations, and
            %   recomputes the LU factorization of the free part of the
            %   coarsest matrix when its mask has or had fixed entries.
            %
            %   OBJ = setupLevelMatrices(OBJ,fixedNodes,computeCoarseLU) also
            %   recomputes the LU factorization when computeCoarseLU is true
            %   (default false).
            if nargin < 3 || isempty(computeCoarseLU)
                computeCoarseLU = false;
            end
            if ~isempty(obj.fixedNodes)
                nnzFixedNodesOld = nnz(obj.fixedNodes{1});
            else
                nnzFixedNodesOld = 0;
            end
            obj.fixedNodes = [];
            obj.fixedNodes{numel(obj.A)} = fixedNodes;
            for curLvl = numel(obj.A) : -1 : 2
                obj.fixedNodes{curLvl-1} = (obj.prolong{curLvl-1}'*obj.fixedNodes{curLvl}) > 0;
            end
            if nnz(obj.fixedNodes{1})>0 || nnzFixedNodesOld>0 || computeCoarseLU
                A0 = obj.A{1};
                L_A=[];U_A=[];P_A=[];Q_A=[];
                if nnz(obj.fixedNodes{1}) ~= numel(obj.fixedNodes{1})
                    [L_A,U_A,P_A,Q_A] = lu(A0(~obj.fixedNodes{1},~obj.fixedNodes{1}));
                end
                obj.L = L_A;
                obj.U = U_A;
                obj.P = P_A;
                obj.Q = Q_A;
            end
        end

        function val = contractionNumber(obj)
            %contractionNumber - Estimate the contraction number
            %   VAL = contractionNumber(OBJ) estimates the spectral radius of
            %   the error propagation matrix I - C*A of one multigrid cycle C
            %   on the finest level with eigs (tolerance 1e-3), assuming a
            %   symmetric cycle.
            curA = obj.A{end};
            curFixedNodes = obj.fixedNodes{end};
            opts.tol = 1e-3;
            opts.issym = true;
            % spectral radius of (I-MG*A)
            val = eigs(@(z)(z-obj.solve(curA(~curFixedNodes,~curFixedNodes)*z,0,1,[],false)),nnz(~curFixedNodes),[],1,'LM',opts);
            val = abs(val);
        end

        function val = dampingFactor(obj)
            %dampingFactor - Estimate the spectral radius of T*A
            %   VAL = dampingFactor(OBJ) estimates the spectral radius of
            %   T*A on the finest level, where the smoother step is
            %   x - T*(A*x - b). The damped smoother converges for
            %   VAL < 2, and VAL close to 1 damps the high frequencies
            %   best, so VAL indicates how to scale the damping.
            curB = zeros(size(obj.A{end},1),1);
            curA = obj.A{end};
            curDamp = obj.damp{end};
            curFixedNodes = obj.fixedNodes{end};
            opts.tol = 1e-2;
            % compute spectral radius of smoother (T*A)
            val = eigs(@(z)(z-obj.smoother(z,curB(~curFixedNodes),curA(~curFixedNodes,~curFixedNodes),curDamp(~curFixedNodes),~curFixedNodes,numel(obj.A))),nnz(~curFixedNodes),[],1,'LM',opts);
            val = abs(val);
        end

        function [x,iter,reduction,contraction] = solve(obj,b,tol,maxIter,x,display)
            %SOLVE - Solve with the multigrid iteration
            %   X = SOLVE(OBJ,B) solves A*X = B for the matrix A of the finest
            %   level restricted to its free entries with multigrid cycles from
            %   the zero initial iterate until the Euclidean norm of the
            %   residual A*X - B is at most 1e-6, with at most 100 cycles. B is
            %   a column vector with nnz(~OBJ.fixedNodes{end}) entries. When a
            %   cycle increases the residual norm, the iteration stops, returns
            %   the previous iterate and prints a message, also for DISPLAY = 0.
            %
            %   X = SOLVE(OBJ,B,TOL) stops when the residual norm is at most
            %   TOL (absolute, not relative to B). TOL must not be empty.
            %
            %   X = SOLVE(OBJ,B,TOL,maxIter) also specifies the maximum number
            %   of cycles (default 100).
            %
            %   X = SOLVE(OBJ,B,...,X0) also specifies the initial iterate, of
            %   the size of B (default zeros).
            %
            %   X = SOLVE(OBJ,B,...,DISPLAY) also selects the output: 0 (no
            %   summary line), 1 (default, one summary line), or 2 (also
            %   computes and prints the contraction number). Empty maxIter,
            %   X0 and DISPLAY select their defaults.
            %
            %   [X,ITER] = SOLVE(...) also returns the number of accepted
            %   cycles.
            %
            %   [X,ITER,REDUCTION] = SOLVE(...) also returns the reduction
            %   factor (R/R0)^(1/(ITER+1)) of the final and the initial
            %   residual norms R and R0.
            %
            %   [X,ITER,REDUCTION,CONTRACTION] = SOLVE(...) also returns the
            %   contraction number of contractionNumber, computed only for
            %   DISPLAY >= 2 (otherwise 0).
            if nargin < 3
                tol = 1e-6;
            end
            if nargin < 4 || isempty(maxIter)
                maxIter = 100;
            end
            if nargin < 5 || isempty(x)
                x = zeros(numel(b),1);
            end
            if nargin < 6 || isempty(display)
                display = true;
            end

            if numel(b)~=nnz(~obj.fixedNodes{end})
                error('fem4m:Multigrid:sizeMismatch','MG: solve input b is not of the same size as the (reduced) system matrix');
            end
            if numel(x)~=nnz(~obj.fixedNodes{end})
                error('fem4m:Multigrid:sizeMismatch','MG: solve input x is not of the same size as the (reduced) system matrix');
            end

            curA = obj.A{end};
            curFixedNodes = obj.fixedNodes{end};
            r = curA(~curFixedNodes,~curFixedNodes)*x-b;
            eNew = sqrt(r'*r);
            iter = 0;
            e0 = eNew;
            while (eNew > tol) && (iter < maxIter)
                iter = iter + 1;
                xOld = x;
                x = Multigrid.MG(x,b,obj.A,obj.smoother,obj.damp,obj.m,obj.prolong,obj.fixedNodes,length(obj.A),obj.L,obj.U,obj.P,obj.Q,obj.w);
                eOld = eNew;
                r = curA(~curFixedNodes,~curFixedNodes)*x-b;
                eNew = sqrt(r'*r);
                if(eNew/eOld > 1)
                    x = xOld;
                    r = curA(~curFixedNodes,~curFixedNodes)*x-b;
                    eNew = sqrt(r'*r);
                    disp(['No convergence, return iterate of step: ',num2str(iter-1)]);
                    iter = iter-1;
                    break;
                end
            end
            reduction = nthroot(eNew/e0,iter+1);
            contraction = 0;
            if display == 1
                disp(['MG: iter = ',int2str(iter),' m = ',num2str(obj.m{end}),' damp = ',num2str(mean(obj.damp{end})),' residual = ',num2str(eNew),' reduction factor = ',num2str(reduction)]);
            elseif display >= 2
                if nnz(~obj.fixedNodes{end}) > 0
                    contraction = obj.contractionNumber;
                end
                disp(['MG: iter = ',int2str(iter),' m = ',num2str(obj.m{end}),' damp = ',num2str(mean(obj.damp{end})),' residual = ',num2str(eNew),' reduction factor = ',num2str(reduction),' contraction number = ',num2str(contraction)]);
            end
        end
    end

    methods(Static)

        function x  = MG(x,b,A,smoother,damp,m,prolong,fixedNodes,level,L,U,P,Q,w)
            %MG - Recursive multigrid cycle
            %   X = MG(X,B,A,SMOOTHER,DAMP,M,PROLONG,fixedNodes,LEVEL,L,U,P,Q,W)
            %   performs one multigrid cycle on level LEVEL for the vectors X
            %   and B on the free entries of that level: M{LEVEL} smoothing
            %   steps, W coarse grid corrections with the restriction
            %   PROLONG{LEVEL-1}' and a recursive cycle on level LEVEL-1, and
            %   M{LEVEL} smoothing steps. A, DAMP, M, PROLONG and fixedNodes
            %   are the cell arrays of the levels (see the properties),
            %   SMOOTHER is the smoother handle, and L, U, P and Q are the LU
            %   factors of the coarsest level. On level 1 the input X is
            %   ignored and the coarse system is solved with the LU factors.
            if level == 1
                x = Q*(U\(L\(P*b)));
            else
                curA = A{level};
                curDamp = damp{level};
                curM = m{level};
                curFixedNodes = fixedNodes{level};
                curP = prolong{level-1};

                for k = 1 : curM
                    x = smoother(x,b,curA(~curFixedNodes,~curFixedNodes),curDamp(~curFixedNodes),~curFixedNodes,level);
                end

                for k = 1 : w
                    defect = curP(~curFixedNodes,~fixedNodes{level-1})'*(b - curA(~curFixedNodes,~curFixedNodes)*x);

                    x0 = zeros(size(defect));
                    correction = Multigrid.MG(x0,defect,A,smoother,damp,m,prolong,fixedNodes,level-1,L,U,P,Q,w);

                    x = x + curP(~curFixedNodes,~fixedNodes{level-1})*correction;
                end

                for k = 1 : curM
                    x = smoother(x,b,curA(~curFixedNodes,~curFixedNodes),curDamp(~curFixedNodes),~curFixedNodes,level);
                end

            end
        end

        function val = Richardson(x,b,A,damp,~,~)
            %RICHARDSON - Richardson smoother
            %   VAL = RICHARDSON(X,B,A,DAMP) applies one damped Richardson step
            %   VAL = X - DAMP.*(A*X - B).
            %
            %   VAL = RICHARDSON(X,B,A,DAMP,FREE,LEVEL) is the call of the
            %   smoother interface (see setupSmoother); FREE and LEVEL are
            %   ignored.
            val = x - damp.*(A*x-b);
        end

        function val = Jacobi(x,b,A,damp,~,~)
            %JACOBI - Jacobi smoother
            %   VAL = JACOBI(X,B,A,DAMP) applies one damped Jacobi step
            %   VAL = X - DAMP./diag(A).*(A*X - B).
            %
            %   VAL = JACOBI(X,B,A,DAMP,FREE,LEVEL) is the call of the smoother
            %   interface (see setupSmoother); FREE and LEVEL are ignored.
            val = x - damp.*1./(diag(A)).*(A*x-b);
        end

        function val = Poisson(x,b,A,damp,~,level,PoissonMG)
            %POISSON - Cycle of a second multigrid solver as smoother
            %   VAL = POISSON(X,B,A,DAMP,FREE,LEVEL,PoissonMG) applies the step
            %   VAL = X - DAMP.*C*(A*X - B), where C is one cycle of the
            %   multigrid solver PoissonMG (e.g. for the Laplacian) on level
            %   LEVEL with the zero initial iterate; FREE is ignored. PoissonMG
            %   must have the same levels and fixed entries as the solver that
            %   uses the smoother. As SMOOTHER it needs a wrapper such as
            %   @(x,b,A,damp,free,level)Multigrid.Poisson(x,b,A,damp,free,level,pmg).
            val = x - damp.*Multigrid.MG(zeros(size(x,1),1),A*x-b,PoissonMG.A,PoissonMG.smoother,PoissonMG.damp,...
                PoissonMG.m,PoissonMG.prolong,PoissonMG.fixedNodes,level,...
                PoissonMG.L,PoissonMG.U,PoissonMG.P,PoissonMG.Q,PoissonMG.w);
        end

    end

end
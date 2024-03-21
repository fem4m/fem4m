classdef Schwarz
    %SCHWARZ - Base class of additive and multiplicative Schwarz methods
    %   Schwarz is the abstract base class of the additive and
    %   multiplicative Schwarz preconditioners; subclasses implement
    %   addLevel to define the subspace decomposition. Create one of the
    %   subclasses, add the levels from the coarsest to the finest mesh
    %   with addLevel, and pass @(R)applyAdd(OBJ,R) or
    %   @(R)applyMult(OBJ,R) as the preconditioner to pcg. The
    %   preconditioner acts on vectors of the free degrees of freedom
    %   (FE.freeNodes) of the finest level added.
    %
    %   Level j stores the system matrix A{j} and the prolongation I{j}
    %   from level j-1, both restricted to the free degrees of freedom
    %   (I{1} = 1), and the damped local solve operator B{j} (a sparse
    %   matrix, or the scalar damping factor for SchwarzBPX and on the
    %   first level of SchwarzHB). With pseudoinverse = true, addLevel
    %   of SchwarzBDDC, SchwarzBPS and SchwarzVO inverts the local
    %   matrices with pinv instead of inv; set it before adding the
    %   levels.
    %
    %   SCHWARZ functions:
    %       addLevel              - Add a decomposition level
    %       applyAdd              - Apply the additive Schwarz preconditioner
    %       applyMult             - Apply the multiplicative preconditioner
    %       additiveSchwarz       - Additive Schwarz algorithm
    %       multiplicativeSchwarz - Multiplicative Schwarz algorithm
    %
    %   SCHWARZ properties:
    %       I             - Prolongations from the next coarser level
    %       B             - Damped local solve operators per level
    %       A             - System matrices per level (free dofs)
    %       omega         - Default damping factor of addLevel
    %       pseudoinverse - Use pinv for the local solves (default false)
    %
    %   See also SchwarzBDDC, SchwarzBPS, SchwarzBPX,
    %   SchwarzHB, SchwarzNN, SchwarzVO

    % Copyright (c) 2026 by the fem4m authors

    properties
        I  % Prolongations from the next coarser level
        B  % Damped local solve operators per level
        A  % System matrices per level (free dofs)
        omega  % Default damping factor of addLevel
        pseudoinverse = false  % Use pinv for the local solves (default false)
    end

    methods (Abstract)
        obj = addLevel(obj,fe,feOld,I,A,A4e,omega);
        %addLevel - Add a decomposition level
        %   OBJ = addLevel(OBJ,FE,FEOLD,I,A) adds the level of the
        %   finite element space FE with the system matrix A on all
        %   FE.nrDoF degrees of freedom. FEOLD is the space of the
        %   previous level and I = prolongation(FE,FEOLD) the
        %   prolongation from it on all degrees of freedom; both are
        %   empty for the first (coarsest) level. The level stores A and
        %   I restricted to FE.freeNodes and the local solve operator,
        %   damped with the property omega. This form is enough for the
        %   subclasses that ignore A4E.
        %
        %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E) also passes the local
        %   element matrices A4E, the second output of
        %   BilinearForm.assemble, which SchwarzBDDC and SchwarzNN need;
        %   the other subclasses ignore A4E.
        %
        %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E,OMEGA) uses the damping
        %   factor OMEGA instead; an empty OMEGA takes omega.
    end

    methods

        function val = applyAdd(obj,b)
            %applyAdd - Apply the additive Schwarz preconditioner
            %   VAL = applyAdd(OBJ,R) applies the additive Schwarz
            %   preconditioner to the vector R on the free degrees of
            %   freedom of the finest level: VAL is the sum over all
            %   levels j of P_j*B{j}*P_j'*R, where P_j is the product of
            %   the prolongations I from level j to the finest level.
            val = Schwarz.additiveSchwarz(b,obj.B,obj.I);
        end

        function val = applyMult(obj,b)
            %applyMult - Apply the multiplicative preconditioner
            %   VAL = applyMult(OBJ,R) applies the symmetric
            %   multiplicative Schwarz preconditioner to the vector R on
            %   the free degrees of freedom of the finest level: on each
            %   level a local solve with B, the correction from the
            %   coarser levels (transferred with I), and a second local
            %   solve with the updated residual (see
            %   multiplicativeSchwarz).
            val = Schwarz.multiplicativeSchwarz(b,obj.B,obj.I,obj.A);
        end

    end

    methods(Static)

        function x = additiveSchwarz(b,B,I,j)
            %additiveSchwarz - Additive Schwarz algorithm
            %   X = additiveSchwarz(R,B,I) applies the additive Schwarz
            %   operator of the levels numel(B),...,1 with the local
            %   solve operators B and the prolongations I (cell arrays)
            %   to the vector R of the finest level J = numel(B):
            %   X = B{J}*R + I{J}*additiveSchwarz(I{J}'*R,B,I,J-1).
            %
            %   X = additiveSchwarz(R,B,I,J) uses the levels J,...,1,
            %   with R on level J.
            if nargin < 4
                j = numel(B);
            end
            x = B{j}*b;
            if j > 1
                c = Schwarz.additiveSchwarz(I{j}'*b,B,I,j-1);
                x = x + I{j}*c;
            end
        end

        function x = multiplicativeSchwarz(b,B,I,A,j)
            %multiplicativeSchwarz - Multiplicative Schwarz algorithm
            %   X = multiplicativeSchwarz(R,B,I,A) applies the
            %   multiplicative Schwarz operator of the levels
            %   numel(B),...,1 with the local solve operators B, the
            %   prolongations I and the system matrices A (cell arrays)
            %   to the vector R of the finest level J = numel(B): first
            %   X = B{J}*R, then the coarse correction
            %   X = X + I{J}*multiplicativeSchwarz(I{J}'*(R-A{J}*X),...),
            %   and for symmetry X = X + B{J}*(R-A{J}*X).
            %
            %   X = multiplicativeSchwarz(R,B,I,A,J) starts on level J,
            %   with R on level J.
            if nargin < 5
                j = numel(B);
            end
            % x0 = 0 -> x0 + B*(b-A*x0)=B*b
            x = B{j}*b;
            % coarse grid correction
            if j > 1
                c  = Schwarz.multiplicativeSchwarz(I{j}'*(b-A{j}*x),B,I,A,j-1);
                x = x + I{j}*c;
            end
            % symmetry
            x = x + B{j}*(b - A{j}*x);
        end

    end
end
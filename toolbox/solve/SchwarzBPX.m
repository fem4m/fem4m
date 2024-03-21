classdef SchwarzBPX < Schwarz
    %SchwarzBPX - BPX multilevel preconditioner
    %   OBJ = SchwarzBPX() creates a BPX (Bramble-Pasciak-Xu)
    %   multilevel preconditioner with the damping factor 1/8.
    %
    %   OBJ = SchwarzBPX(OMEGA) uses the damping factor OMEGA.
    %
    %   The local solve operator of each level is OMEGA times the
    %   identity on its free degrees of freedom (B{j} is the scalar
    %   OMEGA), so applyAdd sums the damped residuals of all levels,
    %   transferred with the prolongations. SchwarzBPX inherits
    %   applyAdd, applyMult and the properties of Schwarz.
    %
    %   SchwarzBPX functions:
    %       addLevel - Add a BPX level
    %
    %   See also Schwarz, SchwarzHB, Multigrid

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = SchwarzBPX(omega)
            %SchwarzBPX - Create BPX preconditioner
            if nargin < 1
                obj.omega = 1/8;
            else
                obj.omega = omega;
            end
        end

        function obj = addLevel(obj,fe,feOld,I,A,~,omega)
            %addLevel - Add a BPX level
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A) adds the level of the
            %   finite element space FE with the full system matrix A,
            %   the space FEOLD of the previous level and the prolongation
            %   I = prolongation(FE,FEOLD), and stores the damping factor
            %   omega as its local solve operator B{j}. An empty I sets the
            %   first level, replacing an existing one; FEOLD is then not
            %   used.
            %
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E) accepts and ignores
            %   the local element matrices A4E.
            %
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E,OMEGA) uses the damping
            %   factor OMEGA instead; an empty OMEGA takes omega.
            if nargin < 7 || isempty(omega)
                omega = obj.omega;
            end
            if isempty(I)
                obj.I{1} = 1;
                obj.B{1} = omega;
                obj.A{1} = A(fe.freeNodes,fe.freeNodes);
            else
                obj.I{end+1} = I(fe.freeNodes,feOld.freeNodes);
                obj.B{end+1} = omega;
                obj.A{end+1} = A(fe.freeNodes,fe.freeNodes);
            end
        end

    end
end
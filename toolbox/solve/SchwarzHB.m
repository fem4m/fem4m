classdef SchwarzHB < Schwarz
    %SchwarzHB - Hierarchical basis preconditioner
    %   OBJ = SchwarzHB() creates a hierarchical basis preconditioner
    %   with the damping factor 1/8.
    %
    %   OBJ = SchwarzHB(OMEGA) uses the damping factor OMEGA.
    %
    %   The local solve operator of the first level is OMEGA times the
    %   identity on all its free degrees of freedom (B{1} is the scalar
    %   OMEGA, no coarse solve), that of each finer level OMEGA times
    %   the identity on its new free degrees of freedom, those at nodes
    %   that are no nodes of the previous level. addLevel finds them in
    %   the prolongation, where a degree of freedom at a node of the
    %   previous level has a unit row, so the numbers of the nodes may
    %   change from level to level. This needs a nodal basis on nested
    %   meshes, as for Pk of any degree on meshes refined with
    %   Mesh.refine. applyMult needs OMEGA*lambda_max(A{j}) < 2 on every
    %   level, which the default damping violates from degree 4 on.
    %   SchwarzHB inherits applyAdd, applyMult and the properties of
    %   Schwarz.
    %
    %   SchwarzHB functions:
    %       addLevel - Add a hierarchical basis level
    %
    %   See also Schwarz, SchwarzBPX, Multigrid

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = SchwarzHB(omega)
            %SchwarzHB - Create hierarchical basis preconditioner
            if nargin < 1
                obj.omega = 1/8;
            else
                obj.omega = omega;
            end
        end

        function obj = addLevel(obj,fe,feOld,I,A,~,omega)
            %addLevel - Add a hierarchical basis level
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A) adds the level of the
            %   finite element space FE with the full system matrix A,
            %   the space FEOLD of the previous level and the prolongation
            %   I = prolongation(FE,FEOLD), and stores the damped identity
            %   on the new free degrees of freedom, those whose row of I
            %   is no unit vector, as its local solve operator. If the
            %   number of unit rows of I differs from the number of
            %   degrees of freedom of FEOLD, as for a basis that is not
            %   nodal, it raises the error fem4m:SchwarzHB:notNested. An
            %   empty I sets the first level, replacing an existing one;
            %   FEOLD is then not used.
            %
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E) accepts and ignores
            %   the local element matrices A4E.
            %
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E,OMEGA) uses the damping
            %   factor OMEGA instead; an empty OMEGA takes omega.
            if nargin < 7 || isempty(omega)
                omega = obj.omega;
            end
            if nargin < 4 || isempty(I)
                obj.I{1} = 1;
                obj.B{1} = omega;
                obj.A{1} = A(fe.freeNodes,fe.freeNodes);
            else
                % a dof at a node of the previous level has a unit row
                % in I; the numbers of the dofs may differ (from P3 on)
                tol = 1e-10;
                atOldNode = full(sum(abs(I) > tol,2) == 1 & abs(max(I,[],2) - 1) < tol);
                if nnz(atOldNode) ~= feOld.nrDoF
                    error('fem4m:SchwarzHB:notNested', ...
                        ['SchwarzHB: %d degrees of freedom of the previous level, ' ...
                        'but %d unit rows in the prolongation.'],feOld.nrDoF,nnz(atOldNode));
                end
                obj.I{end+1} = I(fe.freeNodes,feOld.freeNodes);
                n = numel(fe.freeNodes);
                obj.B{end+1} = omega*spdiags(double(~atOldNode(fe.freeNodes)),0,n,n);
                obj.A{end+1} = A(fe.freeNodes,fe.freeNodes);
            end
        end

    end
end
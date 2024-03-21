classdef SchwarzVO < Schwarz
    %SchwarzVO - Vertex patch preconditioner
    %   OBJ = SchwarzVO() creates an overlapping vertex patch
    %   preconditioner with the damping factor 1/2.
    %
    %   OBJ = SchwarzVO(OMEGA) uses the damping factor OMEGA.
    %
    %   On each level the local solve operator is the sum of the
    %   Dirichlet problems of the system matrix A on the vertex patches,
    %   scaled on both sides with the inverse number of patches that
    %   share a degree of freedom, and damped with OMEGA. The local
    %   problems live on the vertex patches without the dofs on their
    %   boundary; for discontinuous elements, which have no dofs on the
    %   sides, they hold all dofs of the patch elements. SchwarzVO
    %   inherits applyAdd, applyMult and the properties of Schwarz.
    %
    %   SchwarzVO functions:
    %       addLevel - Add a vertex patch level
    %
    %   See also Schwarz, SchwarzBPS, SchwarzNN

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = SchwarzVO(omega)
            %SchwarzVO - Create vertex patch preconditioner
            if nargin < 1
                obj.omega = 1/2;
            else
                obj.omega = omega;
            end
        end

        function obj = addLevel(obj,fe,feOld,I,A,~,omega)
            %addLevel - Add a vertex patch level
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A) adds the level of the
            %   finite element space FE with the full system matrix A,
            %   the space FEOLD of the previous level and the prolongation
            %   I = prolongation(FE,FEOLD) (FEOLD and I empty on the first
            %   level), and stores the vertex patch local solve operator
            %   damped with the property omega.
            %
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E) accepts and ignores
            %   the local element matrices A4E.
            %
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E,OMEGA) uses the damping
            %   factor OMEGA instead; an empty OMEGA takes omega.
            if nargin < 7 || isempty(omega)
                omega = obj.omega;
            end
            obj.A{end+1} = A(fe.freeNodes,fe.freeNodes);
            if ~isempty(feOld) && ~isempty(I)
                obj.I{end+1} = I(fe.freeNodes,feOld.freeNodes);
            else
                obj.I{end+1} = 1;
            end
            mesh = fe.mesh;
            nrNodes = mesh.nrNodes;
            e4nT = mesh.e4n';
            dof4e = fe.dof4e;
            % discontinuous elements have no dofs on the sides
            dof4s = [];
            if isa(mesh,'Mesh3d')
                if isprop(fe,'dof4f')
                    dof4s = fe.dof4f;
                end
                s4e = mesh.f4e;
                n4s = mesh.n4f;
            else
                if isprop(fe,'dof4ed')
                    dof4s = fe.dof4ed;
                end
                s4e = mesh.ed4e;
                n4s = mesh.n4ed;
            end
            %vertex patches
            nrDoF4z = zeros(1,nrNodes);
            dofs = cell(nrNodes,1);
            for curNode = 1 : nrNodes
                % patch elements
                [~,~,omega_z] = find(e4nT(:,curNode));
                % patch boundary sides -> edges/faces
                s4bd_z = s4e(omega_z,:);
                s4bd_z = sort(s4bd_z(:));
                indx = all(n4s(s4bd_z,:)~=curNode,2);
                s4bd_z = s4bd_z(indx);
                % interior dofs of omega_z; all of its dofs for elements
                % without side dofs
                if isempty(dof4s)
                    curDofs = unique(dof4e(omega_z,:));
                else
                    curDofs = setdiff(dof4e(omega_z,:),dof4s(s4bd_z,:));
                end
                dofs{curNode} = curDofs(:)';
                nrDoF4z(curNode) = numel(dofs{curNode});
            end
            nrDoFTypes = unique(nrDoF4z);
            vertexInvA = sparse(size(A,1),size(A,2));
            % solve local Dirichlet problems on patches
            for curNrDoF = nrDoFTypes
                if curNrDoF > 0
                    nodes = find(nrDoF4z==curNrDoF);
                    localInvA = zeros(curNrDoF,curNrDoF,numel(nodes));
                    curdof4z = cell2mat(dofs(nodes));
                    for k = 1 : numel(nodes)
                        localA = A(:,curdof4z(k,:));
                        localInvA(:,:,k) = localA(curdof4z(k,:),:);
                    end
                    if obj.pseudoinverse
                        for k = 1 : size(localInvA,3)
                            localInvA(:,:,k) = pinv(localInvA(:,:,k));
                        end
                    else
                        for k = 1 : size(localInvA,3)
                            localInvA(:,:,k) = inv(localInvA(:,:,k));
                        end
                    end
                    localInvA = permute(localInvA,[3 1 2]);
                    % assemble preconditioner
                    [I,J] = localDoFtoGlobalDoF(curdof4z);
                    vertexInvA = vertexInvA + sparse(I(:),J(:),localInvA(:),size(A,1),size(A,2));
                end
            end
            % scale with the inverse number of patches sharing a dof
            % (Pk in 2d: vertex 1, edge 1/2, element 1/3)
            dofs = cell2mat(cellfun(@transpose,dofs,'UniformOutput',false));
            N = spdiags(1./accumarray(dofs,ones(numel(dofs),1),[fe.nrDoF 1]),0,fe.nrDoF,fe.nrDoF);
            vertexInvA = N*vertexInvA*N;

            obj.B{end+1} = omega*vertexInvA(fe.freeNodes,fe.freeNodes);
        end

    end
end
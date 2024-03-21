classdef SchwarzBPS < Schwarz
    %SchwarzBPS - BPS preconditioner
    %   OBJ = SchwarzBPS() creates a BPS (Bramble-Pasciak-Schatz)
    %   substructuring preconditioner with the damping factor 1/2.
    %
    %   OBJ = SchwarzBPS(OMEGA) uses the damping factor OMEGA.
    %
    %   On each level the local solve operator is the sum of the Jacobi
    %   method diag(A)^(-1) on the vertex dofs and the Dirichlet
    %   problems on the interior degrees of freedom of each edge, each
    %   face (3d) and each element, damped with OMEGA.
    %
    %   The vertex dofs are FE.dof4n and the interior dofs the columns
    %   FE.idof4ed, FE.idof4f and FE.idof4e of dof4ed, dof4f and dof4e
    %   if the element has these properties (as PkVector); otherwise
    %   the vertex dofs must be the global dofs 1:nrNodes and the
    %   interior dofs the last columns of dof4ed, dof4f and dof4e, as
    %   in the numbering of FE.dofManager for elements with one dof per
    %   vertex. The elements must have degrees of freedom on the sides
    %   (edges in 2d, faces in 3d); discontinuous elements raise
    %   fem4m:SchwarzBPS:discontinuousElement. SchwarzBPS inherits
    %   applyAdd, applyMult and the properties of Schwarz.
    %
    %   SchwarzBPS functions:
    %       addLevel - Add a BPS level
    %
    %   See also Schwarz, SchwarzBDDC, SchwarzVO

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = SchwarzBPS(omega)
            %SchwarzBPS - Create BPS preconditioner
            if nargin < 1
                obj.omega = 1/2;
            else
                obj.omega = omega;
            end
        end

        function obj = addLevel(obj,fe,feOld,I,A,~,omega)
            %addLevel - Add a BPS level
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A) adds the level of the
            %   finite element space FE with the full system matrix A,
            %   the space FEOLD of the previous level and the prolongation
            %   I = prolongation(FE,FEOLD) (FEOLD and I empty on the first
            %   level), and stores the BPS local solve operator damped with
            %   the property omega.
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
            % the subdomains need dofs on the sides, shared by neighbouring
            % elements, which discontinuous elements do not have
            if isa(mesh,'Mesh3d')
                sideDoFs = 'dof4f';
            else
                sideDoFs = 'dof4ed';
            end
            if ~isprop(fe,sideDoFs) || isempty(fe.(sideDoFs))
                error('fem4m:SchwarzBPS:discontinuousElement', ...
                    'SchwarzBPS needs degrees of freedom on the sides; %s has none.',class(fe));
            end
            dof4e = fe.dof4e;
            dof4ed = fe.dof4ed;
            n4ed = mesh.n4ed;
            nrEdges = mesh.nrEdges;
            if isa(mesh,'Mesh3d')
                dof4f = fe.dof4f;
                dof4s = dof4f;
                s4e = mesh.f4e;
                nrFaces = mesh.nrFaces;
            else
                dof4s = fe.dof4ed;
                s4e = mesh.ed4e;
            end
            if isprop(fe,'dof4n')
                dof4n = fe.dof4n;
            else
                dof4n = 1:mesh.nrNodes;
            end
            if isprop(fe,'idof4e')
                nrInnerDoF = fe.nrBasis4e;
                idof4e = fe.idof4e;
            else
                nrInnerDoF = numel(setdiff(dof4e(1,:),dof4s(s4e(1,:),:)));
                idof4e = size(dof4e,2)-nrInnerDoF+1:size(dof4e,2);
            end
            if isprop(fe,'idof4ed')
                nrEdgeDoF = fe.nrBasis4ed;
                idof4ed = fe.idof4ed;
            else
                nrEdgeDoF = numel(setdiff(dof4ed(1,:),n4ed(1,:)));
                idof4ed = size(dof4ed,2)-nrEdgeDoF+1:size(dof4ed,2);
            end
            if isa(mesh,'Mesh3d')
                if isprop(fe,'idof4f')
                    nrFaceDoF = fe.nrBasis4f;
                    idof4f = fe.idof4f;
                else
                    nrFaceDoF = numel(setdiff(dof4f(1,:),dof4ed(mesh.ed4f(1,:),:)));
                    idof4f = size(dof4f,2)-nrFaceDoF+1:size(dof4f,2);
                end
            end

            % solve Dirichlet on element
            elementInvA = sparse(size(A,1),size(A,2));
            if nrInnerDoF > 0
                localInvA = zeros(nrInnerDoF,nrInnerDoF,mesh.nrElems);
                dof4inner = dof4e(:,idof4e);
                for k = 1 : mesh.nrElems
                    localA = A(:,dof4inner(k,:));
                    localInvA(:,:,k) = localA(dof4inner(k,:),:);
                end
                if obj.pseudoinverse
                    for k = 1 : mesh.nrElems
                        localInvA(:,:,k) = pinv(localInvA(:,:,k));
                    end
                else
                    for k = 1 : mesh.nrElems
                        localInvA(:,:,k) = inv(localInvA(:,:,k));
                    end
                end
                localInvA = permute(localInvA,[3 1 2]);
                % assemble preconditioner
                [I,J] = localDoFtoGlobalDoF(dof4inner);
                elementInvA = elementInvA + sparse(I(:),J(:),localInvA(:),size(A,1),size(A,2));
            end

            % solve Dirichlet on edges
            edgeInvA = sparse(size(A,1),size(A,2));
            if nrEdgeDoF > 0
                localInvA = zeros(nrEdgeDoF,nrEdgeDoF,nrEdges);
                dof4edInner = dof4ed(:,idof4ed);
                for k = 1 : nrEdges
                    localA = A(:,dof4edInner(k,:));
                    localInvA(:,:,k) = localA(dof4edInner(k,:),:);
                end
                if obj.pseudoinverse
                    for k = 1 :nrEdges
                        localInvA(:,:,k) = pinv(localInvA(:,:,k));
                    end
                else
                    for k = 1 :nrEdges
                        localInvA(:,:,k) = inv(localInvA(:,:,k));
                    end
                end
                localInvA = permute(localInvA,[3 1 2]);
                % assemble preconditioner
                [I,J] = localDoFtoGlobalDoF(dof4edInner);
                edgeInvA = edgeInvA + sparse(I(:),J(:),localInvA(:),size(A,1),size(A,2));
            end

            % solve Dirichlet on faces
            faceInvA = sparse(size(A,1),size(A,2));
            if isa(mesh,'Mesh3d')
                if nrFaceDoF > 0
                    localInvA = zeros(nrFaceDoF,nrFaceDoF,nrFaces);
                    dof4fInner = dof4f(:,idof4f);
                    for k = 1 : nrFaces
                        localA = A(:,dof4fInner(k,:));
                        localInvA(:,:,k) = localA(dof4fInner(k,:),:);
                    end
                    if obj.pseudoinverse
                        for k = 1 : nrFaces
                            localInvA(:,:,k) = pinv(localInvA(:,:,k));
                        end
                    else
                        for k = 1 : nrFaces
                            localInvA(:,:,k) = inv(localInvA(:,:,k));
                        end
                    end
                    localInvA = permute(localInvA,[3 1 2]);
                    % assemble preconditioner
                    [I,J] = localDoFtoGlobalDoF(dof4fInner);
                    faceInvA = faceInvA + sparse(I(:),J(:),localInvA(:),size(A,1),size(A,2));
                end
            end

            % vertices diag(A)^(-1)
            vertexInvA = sparse(dof4n,dof4n,1./diag(A(dof4n,dof4n)),size(A,1),size(A,2));

            %preconditioner
            B = vertexInvA + edgeInvA + faceInvA + elementInvA ;
            obj.B{end+1} = omega*B(fe.freeNodes,fe.freeNodes);

        end

    end
end
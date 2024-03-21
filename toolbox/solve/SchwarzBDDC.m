classdef SchwarzBDDC < Schwarz
    %SchwarzBDDC - BDDC preconditioner
    %   OBJ = SchwarzBDDC() creates a BDDC (balancing domain
    %   decomposition by constraints) preconditioner with the elements
    %   as subdomains and the damping factor 1/2.
    %
    %   OBJ = SchwarzBDDC(OMEGA) uses the damping factor OMEGA.
    %
    %   On each level the local solve operator is the sum of the
    %   Dirichlet problems on the interior degrees of freedom of each
    %   element, the element problems on the remaining non-primal
    %   degrees of freedom with zero primal ones (scaled on both sides
    %   with the inverse number of elements that share a degree of
    %   freedom), and the Jacobi method diag(A)^(-1) on the primal
    %   degrees of freedom, damped with OMEGA. The primal degrees of
    %   freedom are the vertex dofs and, in 3d with edge dofs, one dof
    %   per edge; elements without vertex dofs use one dof per edge
    %   instead, or one per face in 3d without edge dofs.
    %
    %   addLevel needs the local element matrices A4E. The element must
    %   define nrBasis4ed (and nrBasis4n for the vertex dofs; without
    %   it one dof per vertex is assumed). The global numbering must
    %   start with the primal dofs (all vertex dofs, then one dof per
    %   edge, then one per face), and dof4e must list the primal dofs
    %   of each element first and its interior dofs last, as the
    %   numbering of FE.dofManager without edge orientation does for
    %   elements with one component (Pk3d of degree >= 3 reorients its
    %   edge dofs and does not satisfy this). The elements must have
    %   degrees of freedom on the sides (edges in 2d, faces in 3d);
    %   discontinuous elements raise
    %   fem4m:SchwarzBDDC:discontinuousElement. SchwarzBDDC inherits
    %   applyAdd, applyMult and the properties of Schwarz.
    %
    %   SchwarzBDDC functions:
    %       addLevel - Add a BDDC level
    %
    %   SchwarzBDDC properties:
    %       nrFixNodes       - Number of primal dofs per vertex
    %       nrFixEdgeMoments - Number of primal dofs per edge (0 or 1)
    %       nrFixFaceMoments - Number of primal dofs per face (0 or 1)
    %
    %   See also Schwarz, SchwarzBPS, SchwarzNN

    % Copyright (c) 2026 by the fem4m authors

    properties
        nrFixNodes  % Number of primal dofs per vertex
        nrFixEdgeMoments  % Number of primal dofs per edge (0 or 1)
        nrFixFaceMoments  % Number of primal dofs per face (0 or 1)
    end

    methods

        function obj = SchwarzBDDC(omega)
            %SchwarzBDDC - Create BDDC preconditioner
            if nargin < 1
                obj.omega = 1/2;
            else
                obj.omega = omega;
            end
        end

        function obj = addLevel(obj,fe,feOld,I,A,A4e,omega)
            %addLevel - Add a BDDC level
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E) adds the level of
            %   the finite element space FE with the full system matrix A,
            %   the space FEOLD of the previous level, the prolongation
            %   I = prolongation(FE,FEOLD) (FEOLD and I empty on the first
            %   level) and the local element matrices A4E, the second
            %   output of BilinearForm.assemble. It sets nrFixNodes,
            %   nrFixEdgeMoments and nrFixFaceMoments for FE and stores the
            %   BDDC local solve operator damped with the property omega.
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
                error('fem4m:SchwarzBDDC:discontinuousElement', ...
                    'SchwarzBDDC needs degrees of freedom on the sides; %s has none.',class(fe));
            end
            dof4e = fe.dof4e;
            obj.nrFixNodes = 1;
            obj.nrFixFaceMoments = 0;
            obj.nrFixEdgeMoments = 0;
            if isprop(fe,'nrBasis4n')
                obj.nrFixNodes = fe.nrBasis4n;
            end
            if obj.nrFixNodes == 0 && fe.nrBasis4ed == 0
                obj.nrFixFaceMoments = 1;
            elseif obj.nrFixNodes == 0
                obj.nrFixEdgeMoments = 1;
            elseif isa(mesh,'Mesh3d') && fe.nrBasis4ed>0
                obj.nrFixEdgeMoments = 1;
            end
            if isa(mesh,'Mesh3d')
                dof4s = fe.dof4f;
                s4e = mesh.f4e;
                nrCDoF = 4*obj.nrFixNodes+6*obj.nrFixEdgeMoments+4*obj.nrFixFaceMoments;
                nrFaces = mesh.nrFaces;
            else
                dof4s = fe.dof4ed;
                s4e = mesh.ed4e;
                nrCDoF = 3*obj.nrFixNodes+3*obj.nrFixEdgeMoments;
                nrFaces = 0;
            end
            A4e = permute(A4e,[2 3 1]);

            % solve Dirichlet on element
            elementInvA = sparse(size(A,1),size(A,2));
            nrInnerDoF = numel(setdiff(fe.dof4e(1,:),dof4s(s4e(1,:),:)));
            if nrInnerDoF > 0
                localInvA = zeros(nrInnerDoF,nrInnerDoF,mesh.nrElems);
                dof4inner = dof4e(:,end-nrInnerDoF+1:end);
                if obj.pseudoinverse
                    for k = 1 : mesh.nrElems
                        localInvA(:,:,k) = pinv(A4e(end-nrInnerDoF+1:end,end-nrInnerDoF+1:end,k));
                    end
                else
                    for k = 1 : mesh.nrElems
                        localInvA(:,:,k) = inv(A4e(end-nrInnerDoF+1:end,end-nrInnerDoF+1:end,k));
                    end
                end
                localInvA = permute(localInvA,[3 1 2]);
                % assemble preconditioner
                [I,J] = localDoFtoGlobalDoF(dof4inner);
                elementInvA = elementInvA + sparse(I(:),J(:),localInvA(:),size(A,1),size(A,2));
            end

            % solve the element problems with zero primal dofs (vertex,
            % edge or face dofs)
            zeroHcInvA = sparse(size(A,1),size(A,2));
            nrZeroHcDoF = size(dof4e,2)-nrInnerDoF-nrCDoF;
            if nrZeroHcDoF > 0
                localInvA = zeros(nrZeroHcDoF,nrZeroHcDoF,mesh.nrElems);
                dof4zerocorner = dof4e(:,nrCDoF+(1:nrZeroHcDoF));
                if obj.pseudoinverse
                    for k = 1 : mesh.nrElems
                        localInvA(:,:,k) = pinv(A4e(nrCDoF+(1:nrZeroHcDoF),nrCDoF+(1:nrZeroHcDoF),k));
                    end
                else
                    for k = 1 : mesh.nrElems
                        localInvA(:,:,k) = inv(A4e(nrCDoF+(1:nrZeroHcDoF),nrCDoF+(1:nrZeroHcDoF),k));
                    end
                end
                localInvA = permute(localInvA,[3 1 2]);
                % assemble preconditioner
                [I,J] = localDoFtoGlobalDoF(dof4zerocorner);
                zeroHcInvA = zeroHcInvA + sparse(I(:),J(:),localInvA(:),size(A,1),size(A,2));
                % mean value at dofs
                N = spdiags(1./accumarray(dof4e(:),ones(numel(dof4e),1),[fe.nrDoF 1]),0,fe.nrDoF,fe.nrDoF);
                zeroHcInvA = N*zeroHcInvA*N;
            end

            % S0: remaining skeleton (vertices and/or edges and/or faces)
            % primal dofs diag(A)^(-1)
            indices = 1 : mesh.nrNodes*obj.nrFixNodes+mesh.nrEdges*obj.nrFixEdgeMoments+nrFaces*obj.nrFixFaceMoments;
            CInvA = sparse(indices,indices,1./diag(A(indices,indices)),size(A,1),size(A,2));

            %preconditioner
            B = CInvA + zeroHcInvA + elementInvA;
            obj.B{end+1} = omega*B(fe.freeNodes,fe.freeNodes);

        end

    end
end
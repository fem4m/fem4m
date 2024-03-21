classdef SchwarzNN < Schwarz
    %SchwarzNN - Neumann-Neumann preconditioner
    %   OBJ = SchwarzNN() creates a Neumann-Neumann preconditioner with
    %   the elements as subdomains and the damping factor 1/2.
    %
    %   OBJ = SchwarzNN(OMEGA) uses the damping factor OMEGA.
    %
    %   On each level the local solve operator is the sum of the
    %   Dirichlet problems on the interior degrees of freedom of each
    %   element, the Neumann problems on the remaining degrees of
    %   freedom of each element (regularized with the element mass
    %   matrix divided by the element area, the volume in 3d, and
    %   scaled on both sides with the inverse number of elements that
    %   share a degree of freedom), and the Jacobi method diag(A)^(-1)
    %   on the vertex dofs, damped with OMEGA. The local matrices are
    %   always inverted with inv; the property pseudoinverse is
    %   ignored.
    %
    %   addLevel needs the local element matrices A4E. The vertex dofs
    %   must be the global dofs 1:nrNodes and dof4e must list the
    %   interior dofs of each element last, as in the numbering of
    %   FE.dofManager for elements with one dof per vertex. The elements
    %   must have degrees of freedom on the sides (edges in 2d, faces in
    %   3d); discontinuous elements raise
    %   fem4m:SchwarzNN:discontinuousElement. SchwarzNN inherits
    %   applyAdd, applyMult and the properties of Schwarz.
    %
    %   SchwarzNN functions:
    %       addLevel - Add a Neumann-Neumann level
    %
    %   See also Schwarz, SchwarzBPS, SchwarzBDDC

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = SchwarzNN(omega)
            %SchwarzNN - Create Neumann-Neumann preconditioner
            if nargin < 1
                obj.omega = 1/2;
            else
                obj.omega = omega;
            end
        end

        function obj = addLevel(obj,fe,feOld,I,A,A4e,omega)
            %addLevel - Add a Neumann-Neumann level
            %   OBJ = addLevel(OBJ,FE,FEOLD,I,A,A4E) adds the level of
            %   the finite element space FE with the full system matrix A,
            %   the space FEOLD of the previous level, the prolongation
            %   I = prolongation(FE,FEOLD) (FEOLD and I empty on the first
            %   level) and the local element matrices A4E, the second
            %   output of BilinearForm.assemble, and stores the
            %   Neumann-Neumann local solve operator damped with the
            %   property omega.
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
                error('fem4m:SchwarzNN:discontinuousElement', ...
                    'SchwarzNN needs degrees of freedom on the sides; %s has none.',class(fe));
            end
            dof4e = fe.dof4e;
            if isa(mesh,'Mesh3d')
                dof4s = fe.dof4f;
                s4e = mesh.f4e;
            else
                dof4s = fe.dof4ed;
                s4e = mesh.ed4e;
            end
            A4e = permute(A4e,[2 3 1]);

            % solve Dirichlet on element
            elementInvA = sparse(size(A,1),size(A,2));
            nrInnerDoF = numel(setdiff(dof4e(1,:),dof4s(s4e(1,:),:)));
            if nrInnerDoF > 0
                localInvA = zeros(nrInnerDoF,nrInnerDoF,mesh.nrElems);
                dof4inner = dof4e(:,end-nrInnerDoF+1:end);
                for k = 1 : mesh.nrElems
                    localInvA(:,:,k) = inv(A4e(end-nrInnerDoF+1:end,end-nrInnerDoF+1:end,k));
                end
                localInvA = permute(localInvA,[3 1 2]);
                % assemble preconditioner
                [I,J] = localDoFtoGlobalDoF(dof4inner);
                elementInvA = elementInvA + sparse(I(:),J(:),localInvA(:),size(A,1),size(A,2));
            end

            % solve Neumann on boundary of subdomains
            bdryInvA = sparse(size(A,1),size(A,2));
            nrBdryDoF = size(dof4e,2)-nrInnerDoF;
            localInvA = zeros(nrBdryDoF,nrBdryDoF,mesh.nrElems);
            dof4bdry = dof4e(:,1:nrBdryDoF);
            m = BilinearForm('u*v dx');
            [~,M4e] = m.assemble(fe);
            M4e = permute(mesh.area4e.^(-1).*M4e,[2 3 1]);
            A4e = A4e + M4e;
            for k = 1 :mesh.nrElems
                localInvA(:,:,k) = inv(A4e(1:nrBdryDoF,1:nrBdryDoF,k));
            end
            localInvA = permute(localInvA,[3 1 2]);
            % assemble preconditioner
            [I,J] = localDoFtoGlobalDoF(dof4bdry);
            bdryInvA = bdryInvA + sparse(I(:),J(:),localInvA(:),size(A,1),size(A,2));
            % mean value at dofs
            N = spdiags(1./accumarray(dof4e(:),ones(numel(dof4e),1),[fe.nrDoF 1]),0,fe.nrDoF,fe.nrDoF);
            bdryInvA = N*bdryInvA*N;

            % vertices diag(A)^(-1)
            vertexInvA = sparse(1:mesh.nrNodes,1:mesh.nrNodes,1./diag(A(1:mesh.nrNodes,1:mesh.nrNodes)),size(A,1),size(A,2));

            %preconditioner
            B = vertexInvA + bdryInvA + elementInvA;
            obj.B{end+1} = omega*B(fe.freeNodes,fe.freeNodes);

        end

    end
end
classdef DGDkFace < FEFace
    %DGDkFace - Discontinuous Dubiner face space of total degree k
    %   OBJ = DGDkFace(DEGREE,MESH) creates the discontinuous
    %   face space of total polynomial degree DEGREE >= 0 on the
    %   tetrahedral mesh MESH. Its basis functions are the
    %   Dubiner basis of Jacobi polynomials on the reference
    %   triangle (as in DGDk), replicated on both sides of each
    %   face (zero on the second side of a boundary face).
    %
    %   DGDkFace functions:
    %       getFaceBasis    - Evaluate the face basis functions
    %       dofManager      - Assemble the degrees of freedom
    %       getFaceRefBasis - Evaluate the face reference basis
    %
    %   DGDkFace properties:
    %       mesh         - Mesh on which the element is defined
    %       degree       - Polynomial degree of the element
    %       nrBasisFunc  - Number of basis functions per face
    %       dof4f        - Local-to-global dofs per face
    %       nrDoF        - Total number of degrees of freedom
    %       freeNodes    - Indices of the free degrees of freedom
    %       fixedNodes   - Indices of the fixed degrees of freedom
    %
    %   See also FEFace, DGDkEdge, DGDk, DGPkFace

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh          % mesh on which the element is defined
        degree        % polynomial degree of the element
        nrBasisFunc   % number of basis functions per face
        dof4f         % local-to-global dofs per face
        nrDoF         % total number of degrees of freedom
        freeNodes     % indices of the free degrees of freedom
        fixedNodes    % indices of the fixed degrees of freedom
    end

    methods

        function obj = DGDkFace(degree,mesh)
            %DGDkFace - Create a discontinuous Dubiner face space
            obj.degree = degree;
            obj.nrBasisFunc = 3 + 3*(obj.degree-1) + (obj.degree-2)*(obj.degree-1)/2;
            obj.mesh = mesh;
            obj = obj.dofManager;
        end

        function val = getFaceBasis(obj,~,pts_ref,parts)
            %getFaceBasis - Evaluate the face basis functions
            %   VAL = getFaceBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   Dubiner face basis of degree OBJ.degree at the
            %   reference point PTS_REF (1-by-2) on the faces PARTS and
            %   replicates the result on both sides of each face,
            %   the elements e4f(PARTS,1) and e4f(PARTS,2); the
            %   second side of a boundary face is zero. VAL is
            %   numel(PARTS)-by-nrBasisFunc-by-2. PTS is ignored.
            val = obj.getFaceRefBasis(obj.degree,pts_ref);
            val = reshape((val(:)*ones(1,obj.mesh.nrFaces))',obj.mesh.nrFaces,[]);
            val = ones(obj.mesh.nrFaces,1,2).*val;
            val(obj.mesh.e4f(:,2)==0,:,2)=0;
            val = val(parts,:,:);
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) numbers the degrees of
            %   freedom face by face, fixes the degrees of freedom
            %   on the Dirichlet faces f4Db, and frees all others.
            obj.nrDoF = obj.nrBasisFunc*obj.mesh.nrFaces;
            obj.dof4f = reshape(1:obj.nrDoF,obj.nrBasisFunc,obj.mesh.nrFaces)';
            obj.fixedNodes = unique(obj.dof4f(obj.mesh.f4Db,:));
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
        end

    end

    methods (Static)
        %% Face basis helper functions
        function val = getFaceRefBasis(p,x)
            %getFaceRefBasis - Evaluate the face reference basis
            %   VAL = getFaceRefBasis(P,X) evaluates the Dubiner
            %   basis of total degree P at the reference points X
            %   (N-by-2) on the reference triangle. VAL is
            %   N-by-(P+1)(P+2)/2.
            val = zeros(size(x,1),(p^2+3*p+2)/2);
            y = x(:,2);
            x = x(:,1);
            indx = 1;
            PJ = JacobiPolynomials.poly(2*x./(1-y)-1,p+1,0,0,0)';
            t2 = 1-y;
            for i = 0 : p
                PK = JacobiPolynomials.poly(2*y-1,p+1,2*i+1,0,0)';
                t3 = t2.^i;
                for j = 0 : p - i
                    J = PJ(:,i+1);
                    K = PK(:,j+1);
                    val(:,indx) = J.*t3.*K;
                    if i > 0; val(y==1,indx) = 0; end
                    indx = indx + 1;
                end
            end
        end
    end
end
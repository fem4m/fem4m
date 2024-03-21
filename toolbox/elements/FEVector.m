classdef FEVector < FE
    %FEVector - Vector-valued finite element built from scalar elements
    %   OBJ = FEVector(FE1) creates a vector-valued element wrapping the
    %   scalar element FE1.
    %
    %   OBJ = FEVector(FE1,FE2,...) stacks the scalar elements FE1, FE2,
    %   ... as the components of one vector-valued element on the mesh
    %   of FE1 and numbers their degrees of freedom component after
    %   component. The elements must be either all reference-element
    %   based or all physical-element based (isRefElemBased), else an
    %   error is raised; the degree is the largest degree of FE1, FE2, ...
    %
    %   The evaluation routines take the arguments PTS, PTS_REF and PARTS
    %   of the scalar elements: the reference point PTS_REF (1-by-d), its
    %   images PTS in the elements PARTS, used by physical-element based
    %   components, and the element indices PARTS. The ...Basis
    %   routines return one row per element (all elements for
    %   reference-element based components, which ignore PARTS), one
    %   column per basis function and the components in the further
    %   dimensions; the ...Value routines return one row per element of
    %   PARTS.
    %
    %   FEVector functions:
    %       getBasis          - Evaluate the basis functions
    %       getGradBasis      - Evaluate the gradients of the basis
    %       getSymGradBasis   - Evaluate the symmetric gradient of the basis
    %       getDivBasis       - Evaluate the divergence of the basis
    %       getD2Basis        - Evaluate second derivatives of the basis
    %       getLaplaceValue   - Evaluate the Laplacian of the vector function
    %       getRotBasis       - Evaluate the rot of the 2d basis
    %       getCurlBasis      - Evaluate the curl of the 3d basis
    %       getCurlRotBasis   - Evaluate the curl of the rot of the 2d basis
    %       getCurlCurlBasis  - Evaluate the curl-curl of the 3d basis
    %       getSymGradValue   - Evaluate the symmetric gradient
    %       prolongation      - Prolongation matrix from a coarser space
    %       interpolate       - Interpolate a function into the space
    %       plot              - Plot the vector-valued function or a component
    %       plotDisplacement  - Plot the displaced mesh
    %       plotBending       - Plot the bending deformation
    %       plotStreamLines   - Plot stream lines of the vector function
    %       plotVelocity      - Plot a normalized velocity field as quiver
    %       plus              - Merge two vector-valued elements
    %
    %   FEVector properties:
    %       mesh           - Mesh on which the elements are defined
    %       fe             - Cell array of the scalar finite elements
    %       dof4e          - Local-to-global dofs per element
    %       nrDoF          - Number of used degrees of freedom
    %       allnrDoF       - Total number of degrees of freedom
    %       nonDoFNodes    - Unused degrees of freedom of elementwise degrees
    %       freeNodes      - Indices of the free degrees of freedom
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       nrBasisFunc    - Number of basis functions per element
    %       degree         - Largest polynomial degree of the components
    %       isRefElemBased - True if the basis is on the reference element
    %
    %   See also FE, FEMatrix, FESymMatrix, DGPkVector

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh             % mesh on which the elements are defined
        fe               % cell array of the scalar finite elements
        dof4e            % local-to-global dofs per element
        nrDoF            % number of used degrees of freedom
        allnrDoF         % total number of degrees of freedom
        nonDoFNodes      % unused degrees of freedom of elementwise degrees
        freeNodes        % indices of the free degrees of freedom
        fixedNodes       % indices of the fixed degrees of freedom
        nrBasisFunc      % number of basis functions per element
        degree           % largest polynomial degree of the components
        isRefElemBased   % true if the basis is on the reference element
    end

    methods
        function obj = FEVector(varargin)
            %FEVector - Create a vector-valued finite element object
            obj.fe = varargin(:);
            obj.mesh = obj.fe{1}.mesh;
            obj.degree = obj.fe{1}.degree;
            obj.isRefElemBased = obj.fe{1}.isRefElemBased;
            for curFE = 2 : numel(obj.fe)
                obj.degree = max(obj.degree,obj.fe{curFE}.degree);
                if obj.isRefElemBased ~= obj.fe{curFE}.isRefElemBased
                    error('fem4m:FEVector:mixedBasisTypes','FEVector the given FE must be either all reference element based or all physical element based');
                end
            end
            obj = obj.dofManager;
        end

        function val = getBasis(obj,pts,pts_ref,parts)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the vector
            %   basis functions of all scalar elements in OBJ.fe at the
            %   reference point PTS_REF (1-by-d) in the elements PARTS.
            %   VAL(:,J,K) is component K of the J-th basis function: the
            %   basis functions of OBJ.fe{K} form the K-th block of
            %   columns and vanish in the other components. The columns
            %   of the elements appended by plus follow.
            val = permute(obj.fe{1}.getBasis(pts,pts_ref,parts),[1 2 4 3]);
            for curFE = 2 : size(obj.fe,1)
                val(:,end+(1:obj.fe{curFE}.nrBasisFunc),curFE,:) = obj.fe{curFE}.getBasis(pts,pts_ref,parts);
            end
            for curBFE = 2 : size(obj.fe,2)
                val(:,end+(1:obj.fe{1,curBFE}.nrBasisFunc),1,:) = obj.fe{1,curBFE}.getBasis(pts,pts_ref,parts);
                for curFE = 2 : size(obj.fe,1)
                    val(:,end+(1:obj.fe{curFE,curBFE}.nrBasisFunc),curFE,:) = obj.fe{curFE,curBFE}.getBasis(pts,pts_ref,parts);
                end
            end
        end

        function val = getGradBasis(obj,pts,pts_ref,parts)
            %getGradBasis - Evaluate the gradients of the basis functions
            %   VAL = getGradBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   Jacobian matrices of the vector basis functions at the
            %   reference point PTS_REF (1-by-d) in the elements PARTS,
            %   VAL(:,J,K,L) = d_L v_K for the J-th basis function v. If
            %   the number of components equals the number of
            %   derivatives, the last two dimensions are merged column by
            %   column, e.g. [d_x v_1, d_x v_2, d_y v_1, d_y v_2] in 2d.
            val = permute(obj.fe{1}.getGradBasis(pts,pts_ref,parts),[1 2 4 3]);
            for curFE = 2 : size(obj.fe,1)
                val(:,end+(1:obj.fe{curFE}.nrBasisFunc),curFE,:) = obj.fe{curFE}.getGradBasis(pts,pts_ref,parts);
            end
            for curBFE = 2 : size(obj.fe,2)
                val(:,end+(1:obj.fe{1,curBFE}.nrBasisFunc),1,:) = obj.fe{1,curBFE}.getGradBasis(pts,pts_ref,parts);
                for curFE = 2 : size(obj.fe,1)
                    val(:,end+(1:obj.fe{curFE,curBFE}.nrBasisFunc),curFE,:) = obj.fe{curFE,curBFE}.getGradBasis(pts,pts_ref,parts);
                end
            end
            if size(val,3)==size(val,4)
                % flatten the square Jacobian matrix but not tensors
                val = val(:,:,:);
            end
        end

        function val = getSymGradBasis(obj,pts,pts_ref,parts)
            %getSymGradBasis - Evaluate the symmetric gradient of the basis
            %   VAL = getSymGradBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   symmetric gradients (grad(v)+grad(v)')/2 of the vector
            %   basis functions v at the reference point PTS_REF (1-by-d)
            %   in the elements PARTS. The entries are ordered
            %   (e_11,e_22,e_12) in 2d and (e_11,e_22,e_33,e_12,e_13,e_23)
            %   in 3d, as in FormType, without a factor 2 on the
            %   off-diagonal entries.
            grad = obj.getGradBasis(pts,pts_ref,parts);
            grad = reshape(grad,size(grad,1),size(grad,2),sqrt(size(grad,3)),sqrt(size(grad,3)));
            grad = 0.5*(grad+permute(grad,[1 2 4 3]));
            dim = size(grad,3);
            val = zeros(size(grad,1),size(grad,2),dim*(dim+1)/2);
            for j = 1 : dim
                val(:,:,j) = grad(:,:,j,j);
            end
            k = dim+1;
            for i = 1 : dim-1
                for j = i+1 : dim
                    val(:,:,k) = grad(:,:,i,j);
                    k = k + 1;
                end
            end
        end

        function val = getDivBasis(obj,pts,pts_ref,parts)
            %getDivBasis - Evaluate the divergence of the basis functions
            %   VAL = getDivBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   divergences of the vector basis functions at the reference
            %   point PTS_REF (1-by-d) in the elements PARTS as the trace
            %   of their Jacobian matrices.
            grad = obj.getGradBasis(pts,pts_ref,parts);
            grad = reshape(grad,size(grad,1),size(grad,2),sqrt(size(grad,3)),sqrt(size(grad,3)));
            val = zeros(size(grad,1),size(grad,2));
            for k = 1 : size(obj.fe,1)
                val = val + grad(:,:,k,k);
            end
        end

        function val = getD2Basis(obj,pts,pts_ref,parts)
            %getD2Basis - Evaluate second derivatives of the basis functions
            %   VAL = getD2Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the second
            %   derivatives of the vector basis functions at the reference
            %   point PTS_REF (1-by-d) in the elements PARTS. VAL(:,J,K,M)
            %   is the M-th second derivative of component K of the J-th
            %   basis function, ordered [xx yy xy] in 2d and
            %   [xx yy zz xy xz yz] in 3d.
            val = permute(obj.fe{1}.getD2Basis(pts,pts_ref,parts),[1 2 4 3]);
            for curFE = 2 : size(obj.fe,1)
                val(:,end+(1:obj.fe{curFE}.nrBasisFunc),curFE,:) = obj.fe{curFE}.getD2Basis(pts,pts_ref,parts);
            end
            for curBFE = 2 : size(obj.fe,2)
                val(:,end+(1:obj.fe{1,curBFE}.nrBasisFunc),1,:) = obj.fe{1,curBFE}.getD2Basis(pts,pts_ref,parts);
                for curFE = 2 : size(obj.fe,1)
                    val(:,end+(1:obj.fe{curFE,curBFE}.nrBasisFunc),curFE,:) = obj.fe{curFE,curBFE}.getD2Basis(pts,pts_ref,parts);
                end
            end
        end

        function val = getLaplaceValue(obj,pts,pts_ref,parts,c)
            %getLaplaceValue - Evaluate the Laplacian of the vector function
            %   VAL = getLaplaceValue(OBJ,PTS,PTS_REF,PARTS,C) sums the pure
            %   second derivatives of each component of the vector-valued
            %   finite element function with coefficients C at the
            %   reference point PTS_REF (1-by-d) in the elements PARTS,
            %   one column per component.
            d2 = obj.getD2Value(pts,pts_ref,parts,c); % xx, yy, xy or xx, yy, zz, xy, xz, yz
            dim = round((sqrt(8*size(d2,3)+1)-1)/2); % 3 -> 2d, 6 -> 3d
            val = sum(d2(:,:,1:dim),3);
        end

        function val = getRotBasis(obj,pts,pts_ref,parts)
            %getRotBasis - Evaluate the rot of the basis functions
            %   VAL = getRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   rotations d_x v_2 - d_y v_1 of the 2d vector basis
            %   functions v at the reference point PTS_REF (1-by-2) in the
            %   elements PARTS from their Jacobian matrices. Elements that
            %   are not 2d vector fields raise the error
            %   fem4m:FE:rotNotDefined.
            if ~obj.isVectorField(2)
                obj.notDefined('rot');
            end
            grad = obj.getGradBasis(pts,pts_ref,parts); % v1_x, v2_x, v1_y, v2_y
            val = grad(:,:,2) - grad(:,:,3);
        end

        function val = getCurlBasis(obj,pts,pts_ref,parts)
            %getCurlBasis - Evaluate the curl of the basis functions
            %   VAL = getCurlBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   curls of the 3d vector basis functions at the reference
            %   point PTS_REF (1-by-3) in the elements PARTS from their
            %   Jacobian matrices. Elements that are not 3d vector fields
            %   raise the error fem4m:FE:curlNotDefined.
            if ~obj.isVectorField(3)
                obj.notDefined('curl');
            end
            grad = obj.getGradBasis(pts,pts_ref,parts);
            grad = reshape(grad,size(grad,1),size(grad,2),3,3); % (component, derivative)
            val = cat(3,grad(:,:,3,2) - grad(:,:,2,3), ...
                grad(:,:,1,3) - grad(:,:,3,1), ...
                grad(:,:,2,1) - grad(:,:,1,2));
        end

        function val = getCurlRotBasis(obj,pts,pts_ref,parts)
            %getCurlRotBasis - Evaluate the curl of the rot of the basis
            %   VAL = getCurlRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates
            %   curl(rot(v)) = (d_y rot(v), -d_x rot(v)) of the 2d vector
            %   basis functions v at the reference point PTS_REF (1-by-2)
            %   in the elements PARTS from their second derivatives.
            %   Elements that are not 2d vector fields raise the error
            %   fem4m:FE:curlRotNotDefined.
            if ~obj.isVectorField(2)
                obj.notDefined('curlRot');
            end
            d2 = obj.getD2Basis(pts,pts_ref,parts); % xx, yy, xy
            val = cat(3,d2(:,:,2,3) - d2(:,:,1,2),d2(:,:,1,3) - d2(:,:,2,1));
        end

        function val = getCurlCurlBasis(obj,pts,pts_ref,parts)
            %getCurlCurlBasis - Evaluate the curl-curl of the basis
            %   VAL = getCurlCurlBasis(OBJ,PTS,PTS_REF,PARTS) evaluates
            %   curl(curl(v)) = grad(div(v)) - laplace(v) of the 3d vector
            %   basis functions v at the reference point PTS_REF (1-by-3)
            %   in the elements PARTS from their second derivatives.
            %   Elements that are not 3d vector fields raise the error
            %   fem4m:FE:curlCurlNotDefined.
            if ~obj.isVectorField(3)
                obj.notDefined('curlCurl');
            end
            d2 = obj.getD2Basis(pts,pts_ref,parts); % xx, yy, zz, xy, xz, yz
            val = cat(3,d2(:,:,2,4) + d2(:,:,3,5) - d2(:,:,1,2) - d2(:,:,1,3), ...
                d2(:,:,1,4) + d2(:,:,3,6) - d2(:,:,2,1) - d2(:,:,2,3), ...
                d2(:,:,1,5) + d2(:,:,2,6) - d2(:,:,3,1) - d2(:,:,3,2));
        end

        function val = getSymGradValue(obj,pts,pts_ref,parts,c)
            %getSymGradValue - Evaluate the symmetric gradient of the function
            %   VAL = getSymGradValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the
            %   symmetric gradient of the vector-valued finite element
            %   function with coefficients C at the reference point
            %   PTS_REF (1-by-d) in the elements PARTS. The entries are
            %   ordered (e_11,e_22,e_12) in 2d and
            %   (e_11,e_22,e_33,e_12,e_13,e_23) in 3d, without a factor 2
            %   on the off-diagonal entries.
            grad = obj.getGradValue(pts,pts_ref,parts,c);
            grad = reshape(grad,size(grad,1),sqrt(size(grad,2)),sqrt(size(grad,2)));
            grad = 0.5*(grad+permute(grad,[1 3 2]));
            dim = size(grad,3);
            val = zeros(size(grad,1),dim*(dim+1)/2);
            for j = 1 : dim
                val(:,j) = grad(:,j,j);
            end
            k = dim+1;
            for i = 1 : dim-1
                for j = i+1 : dim
                    val(:,k) = grad(:,i,j);
                    k = k + 1;
                end
            end
        end

        function P = prolongation(obj,feOld)
            %PROLONGATION - Prolongation matrix from a coarser space
            %   P = PROLONGATION(OBJ,FEOLD) returns the block diagonal
            %   matrix of the prolongations of the components of OBJ.
            %   FEOLD is the coarse vector-valued element, a cell array of
            %   its components, or one scalar element for all components.
            %
            %   For a sum of elements (plus), e.g. the Mini element
            %   PkVector(1,MESH) + FEVector(Bk(3,MESH),Bk(3,MESH)) or
            %   FEMatrix + CurlBRotk, FEOLD is the same sum on the coarse
            %   mesh, and P is the local L2 projection onto the whole fine
            %   sum space on every element (prolongationMatrix on the
            %   physical elements). The degrees of freedom that are
            %   neither free nor fixed, the unused and the linearly
            %   dependent ones, get zero rows. P is exact if the sum
            %   spaces are nested, as BDMk + CurlBRotk and BDMk + CurlBk
            %   under uniform refinement, although their enrichments
            %   alone are not.
            if size(obj.fe,2) > 1 || obj.nrBasisFunc > sum(cellfun(@(f)f.nrBasisFunc,obj.fe(:)))
                % a sum: several columns of elements or the CurlBRotk of
                % FEMatrix appended to the basis
                if ~isa(feOld,class(obj)) || ~isequal(size(feOld.fe),size(obj.fe)) ...
                        || feOld.nrBasisFunc ~= obj.nrBasisFunc
                    error('fem4m:FEVector:prolongationSumMismatch', ...
                        'The old space of a sum of elements must be the same sum on the coarse mesh.');
                end
                feSum = obj;
                feSum.isAffineEquivalent = false;
                feSum.nonDoFNodes = setdiff(1:obj.allnrDoF,[obj.freeNodes(:);obj.fixedNodes(:)])';
                P = prolongationMatrix(feSum,feOld);
                return
            end
            if isa(feOld,'FEVector')
                feOld = feOld.fe;
            end
            if ~iscell(feOld)
                feOld = repmat({feOld},1,numel(obj.fe));
            end
            blocks = cell(1,numel(obj.fe));
            for curFE = 1 : numel(obj.fe)
                blocks{curFE} = obj.fe{curFE}.prolongation(feOld{curFE});
            end
            P = blkdiag(blocks{:});
        end

        function x = interpolate(obj,func)
            %INTERPOLATE - Interpolate a function into the finite element space
            %   X = INTERPOLATE(OBJ,FUNC) evaluates FUNC at the points of
            %   the reference degrees of freedom c4refDof of every
            %   component of the vector-valued element and returns the
            %   allnrDoF coefficients X. As in project,
            %   FUNC is a function of the physical points, @(x)f(x), or of
            %   the points, reference points and element indices,
            %   @(pts,pts_ref,parts)..., and returns the values on all
            %   elements, one column per component. Components appended by
            %   plus, e.g. the bubbles of the Mini element, are not
            %   interpolated. Components without c4refDof raise an error;
            %   use project instead.
            x = zeros(obj.allnrDoF,1);
            curDofs = 0;
            parts = 1:obj.mesh.nrElems;
            for curFE = 1 : size(obj.fe,1)
                if isprop(obj.fe{curFE},'c4refDof')
                    for curRefDoF = 1 : size(obj.fe{curFE}.c4refDof,1)
                        pts_ref = obj.fe{curFE}.c4refDof(curRefDoF,:);
                        funcVal = FormEval.evalFunction(obj.mesh.point4e(pts_ref),pts_ref,parts,func,[]);
                        x(obj.dof4e(:, curDofs+curRefDoF)) = funcVal(:,curFE);
                    end
                    curDofs = curDofs + obj.fe{curFE}.nrBasisFunc;
                else
                    error('fem4m:FEVector:interpolationNotSupported','Interpolation not supported for this vector finite element.\n Use a projection instead.');
                end
            end
        end

        function plot(obj,x,component,varargin)
            %PLOT - Plot the vector-valued function or a component
            %   PLOT(OBJ,X) plots the finite element function with
            %   coefficients X by FE.plot, colored by the Euclidean norm,
            %   or by the von Mises stress for 3 or 4 components in 2d and
            %   for 6 components in 3d.
            %
            %   PLOT(OBJ,X,COMPONENT) plots only the component COMPONENT,
            %   the function of the scalar element OBJ.fe{COMPONENT}, by
            %   the plot method of that element.
            %
            %   PLOT(OBJ,X,COMPONENT,ARGS...) passes the further arguments
            %   ARGS, e.g. the refinement level SUB of FE.plot, to the
            %   plot; COMPONENT = [] plots the whole function.
            if nargin < 3 || isempty(component)
                plot@FE(obj,x,varargin{:});
            else
                % coefficients of the component by the local-to-global maps
                feComponent = obj.fe{component};
                cols = sum(cellfun(@(f)f.nrBasisFunc,obj.fe(1:component-1))) + (1:feComponent.nrBasisFunc);
                xComponent = zeros(max([feComponent.nrDoF;double(feComponent.dof4e(:))]),1);
                xComponent(feComponent.dof4e) = x(obj.dof4e(:,cols));
                feComponent.plot(xComponent,varargin{:});
            end
        end

        function plotDisplacement(obj,c,scale,c2,fe2)
            %plotDisplacement - Plot the displaced mesh
            %   plotDisplacement(OBJ,C) plots, as white triangles, the 2d
            %   mesh displaced by the vector-valued finite element function
            %   with coefficients C, evaluated at the element vertices.
            %
            %   plotDisplacement(OBJ,C,SCALE) amplifies the displacement by
            %   the factor SCALE (default 1).
            %
            %   plotDisplacement(OBJ,C,SCALE,C2,FE2) colors the displaced
            %   mesh by the finite element function of FE2 with
            %   coefficients C2 at the vertices: by its value if it is
            %   scalar, by the plane-stress von Mises stress for 3
            %   components [s11 s22 s12] or 4 components
            %   [s11 s21 s12 s22], and by the Euclidean norm otherwise.
            if nargin < 3 || isempty(scale)
                scale = 1;
            end
            c4n = obj.mesh.c4n;
            n4e = obj.mesh.n4e;
            u = scale*[obj.getValue([],[0 0],[],c),...
                obj.getValue([],[1 0],[],c),...
                obj.getValue([],[0 1],[],c)];
            coord = [c4n(n4e(:,1),:),c4n(n4e(:,2),:),c4n(n4e(:,3),:)];
            if nargin < 5
                patch(coord(:,[1 3 5])'+u(:,[1 3 5])',coord(:,[2 4 6])'+u(:,[2 4 6])','w');
            else
                nodes = [0 0 ; 1 0 ; 0 1];
                color = zeros(obj.mesh.nrElems,size(nodes,1));
                for j = 1 : size(nodes,1)
                    value = fe2.getValue([],nodes(j,:),[],c2);
                    if size(value,2) == 1
                        % value
                        color(:,j) = value;
                    elseif size(value,2) == 3
                        % von Mises stress
                        color(:,j) = sqrt(value(:,1).^2+value(:,2).^2-value(:,1).*value(:,2)+3*value(:,3).^2);
                    elseif size(value,2) == 4
                        % von Mises stress
                        color(:,j) = sqrt(value(:,1).^2+value(:,4).^2-value(:,1).*value(:,4)+3*(1/2*(value(:,2)+value(:,3))).^2);
                    else
                        %norm
                        color(:,j) = sqrt(sum(value.^2,[2 3]));
                    end
                end
                coord = [c4n(n4e(:,1),:),c4n(n4e(:,2),:),c4n(n4e(:,3),:)];
                patch(coord(:,[1 3 5])'+u(:,[1 3 5])',coord(:,[2 4 6])'+u(:,[2 4 6])',color');
            end
        end

        function plotBending(obj,c)
            %plotBending - Plot the bending deformation
            %   plotBending(OBJ,C) plots the surface traced by the
            %   vector-valued finite element function with three
            %   components and coefficients C on a 2d mesh: its values at
            %   the vertices of every triangle are the x, y and z
            %   coordinates of a patch, colored by the third component.
            u = [obj.getValue([],[0 0],[],c),...
                obj.getValue([],[1 0],[],c),...
                obj.getValue([],[0 1],[],c)];
            patch(u(:,[1 4 7])',u(:,[2 5 8])',u(:,[3 6 9])',u(:,[3 6 9])')
        end

        function plotStreamLines(obj,c,anz,ax,bx,ay,by)
            %plotStreamLines - Plot stream lines of the vector function
            %   plotStreamLines(OBJ,C) plots stream lines (streamslice) of
            %   the vector-valued finite element function with
            %   coefficients C on a 2d mesh, interpolated (griddata) from
            %   its values at the element vertices to a grid over the
            %   bounding box of the mesh.
            %
            %   plotStreamLines(OBJ,C,ANZ) uses ANZ intervals, i.e. ANZ+1
            %   grid points, per coordinate direction (default 100).
            %
            %   plotStreamLines(OBJ,C,ANZ,AX,BX,AY,BY) restricts the grid and
            %   the axes to the box [AX,BX] x [AY,BY]; if BY is missing or
            %   empty, the bounding box of the mesh is used.
            if nargin < 3 || isempty(anz)
                anz = 100;
            end
            if nargin < 7  || isempty(by)
                ax = min(obj.mesh.c4n(:,1));
                bx = max(obj.mesh.c4n(:,1));
                ay = min(obj.mesh.c4n(:,2));
                by = max(obj.mesh.c4n(:,2));
            end

            coordX = reshape(obj.mesh.c4n(obj.mesh.n4e,1),[],3);
            coordY = reshape(obj.mesh.c4n(obj.mesh.n4e,2),[],3);
            u4e = [obj.getValue([],[0 0],[],c),...
                obj.getValue([],[1 0],[],c),...
                obj.getValue([],[0 1],[],c)];
            uX = u4e(:,1:2:6);
            uY = u4e(:,2:2:6);
            [XI,YI] = meshgrid(ax:(bx-ax)/anz:bx,ay:(by-ay)/anz:by);
            warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
            UI = griddata(coordX(:),coordY(:),uX(:),XI,YI);
            VI = griddata(coordX(:),coordY(:),uY(:),XI,YI);
            warning('on','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
            h = streamslice(XI,YI,UI,VI);
            set(h,'color','k');
            axis([ax bx ay by]);
        end

        function plotVelocity(obj,c,anz,ax,bx,ay,by)
            %plotVelocity - Plot a normalized velocity field as quiver
            %   plotVelocity(OBJ,C) plots the direction field of the
            %   vector-valued finite element function with coefficients C
            %   on a 2d mesh as unit arrows (quiver) on a grid over the
            %   bounding box of the mesh, interpolated (griddata) from its
            %   values at the element vertices.
            %
            %   plotVelocity(OBJ,C,ANZ) uses ANZ intervals, i.e. ANZ+1 grid
            %   points, per coordinate direction (default 20).
            %
            %   plotVelocity(OBJ,C,ANZ,AX,BX,AY,BY) restricts the grid and the
            %   axes to the box [AX,BX] x [AY,BY]; if BY is missing or
            %   empty, the bounding box of the mesh is used.
            if nargin < 3  || isempty(anz)
                anz = 20;
            end
            if nargin < 7  || isempty(by)
                ax = min(obj.mesh.c4n(:,1));
                bx = max(obj.mesh.c4n(:,1));
                ay = min(obj.mesh.c4n(:,2));
                by = max(obj.mesh.c4n(:,2));
            end
            coordX = reshape(obj.mesh.c4n(obj.mesh.n4e,1),[],3);
            coordY = reshape(obj.mesh.c4n(obj.mesh.n4e,2),[],3);
            u4e = [obj.getValue([],[0 0],[],c),...
                obj.getValue([],[1 0],[],c),...
                obj.getValue([],[0 1],[],c)];
            uX = u4e(:,1:2:6);
            uY = u4e(:,2:2:6);
            [XI,YI] = meshgrid(ax:(bx-ax)/anz:bx,ay:(by-ay)/anz:by);
            warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
            UI = griddata(coordX(:),coordY(:),uX(:),XI,YI);
            VI = griddata(coordX(:),coordY(:),uY(:),XI,YI);
            warning('on','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
            norm = sqrt(UI.^2+VI.^2);
            UI=UI./norm;
            VI=VI./norm;
            quiver(XI,YI,UI,VI);
            axis([ax bx ay by]);
        end

        function obj = plus(obj1,obj2)
            %PLUS - Merge two vector-valued finite element objects
            %   OBJ3 = PLUS(OBJ1,OBJ2), or OBJ3 = OBJ1 + OBJ2, appends the
            %   components of OBJ2 as further columns of OBJ3.fe and shifts
            %   the degrees of freedom of OBJ2 by OBJ1.allnrDoF. Both must
            %   have the same number of components, and those of OBJ2 must
            %   be DGDk elements (or subclasses such as CurlBk), else an
            %   error is raised. The degrees of freedom of OBJ2 of degree
            %   at most OBJ1.degree (OBJ1.degree-1 if OBJ1 is built from
            %   RTk), which are linearly dependent on OBJ1, are removed
            %   from the free set.

            if size(obj1.fe,1) ~= size(obj2.fe,1)
                error('fem4m:FEVector:dimensionMismatch','Dimension mismatch of the two FE vectors, addition not possible!')
            end

            isDGDkVector = true;
            for k = 1 : size(obj2.fe,1)
                isDGDkVector = isDGDkVector & isa(obj2.fe{k},'DGDk');
            end
            if ~isDGDkVector
                error('fem4m:FEVector:invalidSummand','For addition of two FE vectors, the second input must be a ''DGDk'' vector!');
            end

            obj = obj1;
            obj.degree = max(obj1.degree,obj2.degree);
            obj.dof4e = [obj1.dof4e,obj1.allnrDoF+obj2.dof4e];
            obj.nrDoF = obj1.nrDoF+obj2.nrDoF;
            obj.allnrDoF = obj1.allnrDoF+obj2.allnrDoF;
            obj.freeNodes = [obj1.freeNodes(:);obj1.allnrDoF+obj2.freeNodes(:)];
            obj.fixedNodes = [obj1.fixedNodes(:);obj1.allnrDoF+obj2.fixedNodes(:)];
            obj.nonDoFNodes = [obj1.nonDoFNodes(:);obj1.allnrDoF+obj2.nonDoFNodes(:)];
            degree4DoF = zeros(1,obj2.nrBasisFunc);
            indx = 0;
            for k = 1 : size(obj2.fe,1)
                degree4DoF(indx+(1:obj2.fe{k}.nrBasisFunc)) = obj2.fe{k}.degree4DoF;
                indx = indx + obj2.fe{k}.nrBasisFunc;
            end
            if isa(obj1.fe{1},'RTk')
                linearDependentDoFs = obj1.allnrDoF+obj2.dof4e(:,degree4DoF<=obj1.degree-1);
            else
                linearDependentDoFs = obj1.allnrDoF+obj2.dof4e(:,degree4DoF<=obj1.degree);
            end
            obj.freeNodes = setdiff(obj.freeNodes,linearDependentDoFs);
            obj.nrBasisFunc = obj1.nrBasisFunc+obj2.nrBasisFunc;
            indx = size(obj1.fe,2);
            for i = 1 : size(obj1.fe,1)
                for k = 1 : size(obj2.fe,2)
                    obj.fe{i,indx+k} = obj2.fe{i,k};
                end
            end

        end

    end

    methods (Access = protected)

        function obj = dofManager(obj)
            %dofManager - Assemble the combined degrees of freedom
            %   OBJ = dofManager(OBJ) concatenates the local-to-global maps,
            %   the free, fixed and unused degrees of freedom, and the basis
            %   counts of all scalar elements in OBJ.fe. Every element is
            %   shifted by the number of degrees of freedom of the elements
            %   before it, for elementwise degrees by allnrDoF, which
            %   includes the unused ones (nonDoFNodes).
            obj.dof4e = zeros(obj.mesh.nrElems,0);
            obj.freeNodes = zeros(0,1);
            obj.fixedNodes = zeros(0,1);
            obj.nonDoFNodes = zeros(0,1);
            obj.nrDoF = 0;
            obj.allnrDoF = 0;
            obj.nrBasisFunc = 0;
            for curFE = 1 : numel(obj.fe)
                scalarFE = obj.fe{curFE};
                obj.dof4e = [obj.dof4e,obj.allnrDoF+scalarFE.dof4e];
                obj.freeNodes = [obj.freeNodes;obj.allnrDoF+scalarFE.freeNodes(:)];
                obj.fixedNodes = [obj.fixedNodes;obj.allnrDoF+scalarFE.fixedNodes(:)];
                if isprop(scalarFE,'allnrDoF')
                    obj.nonDoFNodes = [obj.nonDoFNodes;obj.allnrDoF+scalarFE.nonDoFNodes(:)];
                    obj.allnrDoF = obj.allnrDoF + scalarFE.allnrDoF;
                else
                    obj.allnrDoF = obj.allnrDoF + scalarFE.nrDoF;
                end
                obj.nrDoF = obj.nrDoF + scalarFE.nrDoF;
                obj.nrBasisFunc = obj.nrBasisFunc + scalarFE.nrBasisFunc;
            end
        end

    end

    methods (Access = protected)

        function tf = isMatrixValued(~)
            %isMatrixValued - True for matrix-valued subclasses
            %   TF = isMatrixValued(OBJ) is false for vector fields;
            %   FEMatrix, FESymMatrix and FESkewSymMatrix override it.
            tf = false;
        end

        function tf = isVectorField(obj,dim)
            %isVectorField - True for a vector field of dimension DIM
            %   TF = isVectorField(OBJ,DIM) is true if OBJ is a vector
            %   field with DIM components on a DIM-dimensional mesh.
            tf = ~obj.isMatrixValued() && size(obj.fe,1) == dim ...
                && isa(obj.mesh,'Mesh3d') == (dim == 3);
        end

    end
end
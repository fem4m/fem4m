classdef FE
    %FE - Abstract base class of the element-based finite element objects
    %   FE defines the common interface of the finite element spaces on
    %   the elements of a mesh (FEEdge and FEFace are their counterparts
    %   on the edges and faces; the face space DGPkFace derives from Pk3d
    %   instead), i.e. the mesh and degree-of-freedom information together
    %   with the evaluation routines of basis functions and their
    %   derivatives. Every element constructor takes (DEGREE,MESH), so
    %   that for a uniform degree feval(class(OBJ),OBJ.degree,MESH) builds
    %   the same space on another mesh. Only FEVector, FEMatrix,
    %   FESymMatrix and FESkewSymMatrix themselves are built from
    %   component spaces instead.
    %
    %   The evaluation routines take the physical points PTS (one row per
    %   element of PARTS), the reference point PTS_REF (1-by-d, where d is
    %   the dimension of the reference element: 1-by-2 on triangles, also
    %   on surface meshes, 1-by-3 on tetrahedra) and the element indices
    %   PARTS, as integrate passes them. Elements whose basis lives on
    %   the reference element (isRefElemBased) evaluate at PTS_REF, the
    %   others (e.g. PDGk) at PTS (PDGk3d at the images of PTS_REF if PTS
    %   is empty); on polygonal meshes PTS_REF is empty.
    %   The ...Value routines return one row per element of PARTS, and an
    %   empty PARTS with a nonempty PTS_REF selects all elements; the
    %   ...Basis routines of reference-element based elements may return
    %   all elements (see FormBasis.elementRows). getValue, getGradValue
    %   and getD2Value accept several columns of C (one function per
    %   column); the other ...Value routines use a single coefficient
    %   vector C.
    %
    %   FE functions:
    %       getValue           - Evaluate the finite element function
    %       getGradValue       - Evaluate the gradient of the function
    %       getDivValue        - Evaluate the divergence of the function
    %       getCurlBasis       - Evaluate the curl of the basis functions
    %       getRotBasis        - Evaluate the rot of the basis functions
    %       getCurlCurlBasis   - Evaluate the curl-curl of the basis functions
    %       getCurlRotBasis    - Evaluate the curl of the rot of the basis
    %       getRotRotBasis     - Evaluate the rot of the rot of the basis
    %       getCurlValue       - Evaluate the curl of the function
    %       getRotValue        - Evaluate the rot of the function
    %       getD2Value         - Evaluate the Hessian of the function
    %       getCurlCurlValue   - Evaluate the curl-curl of the function
    %       getCurlRotValue    - Evaluate the curl of the rot of the function
    %       getRotRotValue     - Evaluate the rot of the rot of the function
    %       getDivDivValue     - Evaluate the double divergence of the function
    %       getD3Value         - Evaluate third derivatives of the function
    %       getLaplaceValue    - Evaluate the Laplacian of the function
    %       getBiharmonicValue - Evaluate the biharmonic of the function
    %       getTraceValue      - Evaluate the trace of the function
    %       interpolate        - Interpolate a function into the space
    %       project            - L2-project a function into the space
    %       prolongation       - Prolongation matrix from a coarser space
    %       prolongate         - Prolongate coefficients of a coarser space
    %       plot               - Plot the finite element function
    %
    %   FE properties:
    %       mesh               - Mesh on which the element is defined
    %       dof4e              - Local-to-global dofs per element
    %       nrBasisFunc        - Number of basis functions per element
    %       nrDoF              - Number of used degrees of freedom
    %       isRefElemBased     - True if the basis is on the reference element
    %       isAffineEquivalent - True if the basis is affine equivalent
    %
    %   See also FEVector, FEMatrix, FEEdge, FEFace

    % Copyright (c) 2026 by the fem4m authors

    % members defined by the subclasses (optional ones are checked
    % with isprop or ismethod)
    %#ok<*MCNPN,*MCNPR>

    properties(Abstract)
        mesh             % mesh on which the element is defined
        dof4e            % local-to-global dofs per element
        nrBasisFunc      % number of basis functions per element
        nrDoF            % number of used degrees of freedom
        isRefElemBased   % true if the basis is on the reference element
    end

    % isAffineEquivalent is true if the basis on an element depends on the
    % mesh only through the affine map from the reference element and the
    % orientation of the element (vertex ids and, in 2d, boundary edges);
    % elements with geometry dependent basis coefficients set it to false
    % in their constructor
    properties
        isAffineEquivalent = true % true if the basis is affine equivalent
    end

    methods(Access=public)

        function val = getValue(obj,pts,pts_ref,parts,c)
            %getValue - Evaluate the finite element function
            %   VAL = getValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the finite
            %   element function with coefficients C in the elements PARTS at
            %   the reference point PTS_REF (1-by-d), or at the points PTS for
            %   elements that are not reference-element based (see FE). VAL has
            %   one row per element of PARTS.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            curDoF4e = obj.dof4e(parts,:);
            u4e = reshape(c(curDoF4e(:),:),[size(curDoF4e,1) obj.nrBasisFunc size(c,2)]);
            evalBasis = FormBasis.elementRows(obj.getBasis(pts,pts_ref,parts),parts,obj);
            val = permute(sum(u4e.*evalBasis,2),[1 3:numel(size(evalBasis)) 2]);
        end

        function val = getGradValue(obj,pts,pts_ref,parts,c)
            %getGradValue - Evaluate the gradient of the finite element function
            %   VAL = getGradValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the
            %   gradient of the finite element function with coefficients C in
            %   the elements PARTS at the reference point PTS_REF (1-by-d), or at
            %   the points PTS for elements that are not reference-element based.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            curDoF4e = obj.dof4e(parts,:);
            u4e = reshape(c(curDoF4e(:),:),[size(curDoF4e,1) obj.nrBasisFunc size(c,2)]);
            grad = FormBasis.elementRows(obj.getGradBasis(pts,pts_ref,parts),parts,obj);
            val = permute(sum(u4e.*grad,2),[1 3:numel(size(grad)) 2]);
        end

        function val = getDivValue(obj,pts,pts_ref,parts,c)
            %getDivValue - Evaluate the divergence of the finite element function
            %   VAL = getDivValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the
            %   divergence of the finite element function with coefficients C
            %   in the elements PARTS at the reference point PTS_REF (1-by-d), or
            %   at the points PTS for elements that are not reference-element
            %   based, for elements that provide getDivBasis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            curDoF4e = obj.dof4e(parts,:);
            coeff4e = reshape(c(curDoF4e),size(curDoF4e));
            basisDiv4e = FormBasis.elementRows(obj.getDivBasis(pts,pts_ref,parts),parts,obj);
            val = permute(sum(basisDiv4e.*coeff4e,2),[1 3 2]);
        end

        function val = getCurlBasis(obj,pts,pts_ref,parts)
            %getCurlBasis - Evaluate the curl of the basis functions
            %   VAL = getCurlBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   vector curl (d_y phi, -d_x phi) of the basis functions phi
            %   of a scalar element on a 2d mesh from their gradients.
            %   Vector-valued elements in 3d override it; all other
            %   elements raise the error fem4m:FE:curlNotDefined.
            if ~isa(obj.mesh,'Mesh3d') && ismethod(obj,'getGradBasis')
                grad = obj.getGradBasis(pts,pts_ref,parts);
                if ndims(grad) == 3 && size(grad,3) == 2
                    val = cat(3,grad(:,:,2),-grad(:,:,1));
                    return
                end
            end
            val = [];
            obj.notDefined('curl');
        end

        function val = getRotBasis(obj,~,~,~)
            %getRotBasis - Evaluate the rot of the basis functions
            %   VAL = getRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   rotation d_x v_2 - d_y v_1 of 2d vector-valued basis
            %   functions (row-wise for matrix-valued ones). Elements that
            %   do not override it raise the error fem4m:FE:rotNotDefined.
            val = [];
            obj.notDefined('rot');
        end

        function val = getCurlCurlBasis(obj,pts,pts_ref,parts)
            %getCurlCurlBasis - Evaluate the curl-curl of the basis functions
            %   VAL = getCurlCurlBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   Airy matrices curl(curl(phi)) = [phi_yy -phi_xy; -phi_xy
            %   phi_xx] of the basis functions phi of a scalar element on
            %   a 2d mesh from their Hessians. Vector-valued elements in 3d
            %   override it; all other elements raise the error
            %   fem4m:FE:curlCurlNotDefined.
            if ~isa(obj.mesh,'Mesh3d') && ismethod(obj,'getD2Basis')
                d2 = obj.getD2Basis(pts,pts_ref,parts); % xx, yy, xy
                if ndims(d2) == 3 && size(d2,3) == 3
                    val = d2(:,:,[2 3 3 1]);
                    val(:,:,[2 3]) = -val(:,:,[2 3]);
                    return
                end
            end
            val = [];
            obj.notDefined('curlCurl');
        end

        function val = getCurlRotBasis(obj,~,~,~)
            %getCurlRotBasis - Evaluate the curl of the rot of the basis
            %   VAL = getCurlRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates
            %   curl(rot(v)) = (d_y rot(v), -d_x rot(v)) of 2d
            %   vector-valued basis functions v. Elements that do not
            %   override it raise the error fem4m:FE:curlRotNotDefined.
            val = [];
            obj.notDefined('curlRot');
        end

        function val = getRotRotBasis(obj,~,~,~)
            %getRotRotBasis - Evaluate the rot of the rot of the basis
            %   VAL = getRotRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   scalar rot(rot(sigma)) of 2d symmetric matrix-valued basis
            %   functions sigma. Elements that do not override it raise
            %   the error fem4m:FE:rotRotNotDefined.
            val = [];
            obj.notDefined('rotRot');
        end

        function val = getCurlValue(obj,pts,pts_ref,parts,c)
            %getCurlValue - Evaluate the curl of the finite element function
            %   VAL = getCurlValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the curl
            %   of the finite element function with coefficients C in the
            %   elements PARTS at the reference point PTS_REF (1-by-d), or at the
            %   points PTS for elements that are not reference-element based,
            %   from getCurlBasis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            val = obj.combineBasis(obj.getCurlBasis(pts,pts_ref,parts),parts,c);
        end

        function val = getRotValue(obj,pts,pts_ref,parts,c)
            %getRotValue - Evaluate the rot of the finite element function
            %   VAL = getRotValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the
            %   rotation of the finite element function with coefficients C
            %   in the elements PARTS at the reference point PTS_REF (1-by-d), or
            %   at the points PTS for elements that are not reference-element
            %   based, from getRotBasis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            val = obj.combineBasis(obj.getRotBasis(pts,pts_ref,parts),parts,c);
        end

        function val = getD2Value(obj,pts,pts_ref,parts,c)
            %getD2Value - Evaluate the Hessian of the finite element function
            %   VAL = getD2Value(OBJ,PTS,PTS_REF,PARTS,C) evaluates the second
            %   derivatives of the finite element function with coefficients C
            %   in the elements PARTS at the reference point PTS_REF (1-by-d), or
            %   at the points PTS for elements that are not reference-element
            %   based, in the component order of getD2Basis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            curDoF4e = obj.dof4e(parts,:);
            u4e = reshape(c(curDoF4e(:),:),[size(curDoF4e,1) obj.nrBasisFunc size(c,2)]);
            d2 = FormBasis.elementRows(obj.getD2Basis(pts,pts_ref,parts),parts,obj);
            val = permute(sum(u4e.*d2,2),[1 3:numel(size(d2)) 2]);
        end

        function val = getCurlCurlValue(obj,pts,pts_ref,parts,c)
            %getCurlCurlValue - Evaluate the curl-curl of the finite element function
            %   VAL = getCurlCurlValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the
            %   curl-curl of the finite element function with coefficients C
            %   in the elements PARTS at the reference point PTS_REF (1-by-d), or
            %   at the points PTS for elements that are not reference-element
            %   based, from getCurlCurlBasis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            val = obj.combineBasis(obj.getCurlCurlBasis(pts,pts_ref,parts),parts,c);
        end

        function val = getCurlRotValue(obj,pts,pts_ref,parts,c)
            %getCurlRotValue - Evaluate the curl of the rot of the function
            %   VAL = getCurlRotValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates
            %   curl(rot(v)) of the 2d vector-valued finite element function v
            %   with coefficients C in the elements PARTS at the reference point
            %   PTS_REF (1-by-d), or at the points PTS for elements that are not
            %   reference-element based, from getCurlRotBasis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            val = obj.combineBasis(obj.getCurlRotBasis(pts,pts_ref,parts),parts,c);
        end

        function val = getRotRotValue(obj,pts,pts_ref,parts,c)
            %getRotRotValue - Evaluate the rot of the rot of the function
            %   VAL = getRotRotValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates
            %   rot(rot(sigma)) of the 2d symmetric matrix-valued finite element
            %   function sigma with coefficients C in the elements PARTS at the
            %   reference point PTS_REF (1-by-d), or at the points PTS for
            %   elements that are not reference-element based, from
            %   getRotRotBasis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            val = obj.combineBasis(obj.getRotRotBasis(pts,pts_ref,parts),parts,c);
        end

        function val = getDivDivValue(obj,pts,pts_ref,parts,c)
            %getDivDivValue - Evaluate the double divergence of the function
            %   VAL = getDivDivValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates
            %   div(div(sigma)) of the symmetric matrix-valued finite
            %   element function sigma with coefficients C in the elements PARTS
            %   at the reference point PTS_REF (1-by-d), or at the points PTS for
            %   elements that are not reference-element based, for elements that
            %   provide getDivDivBasis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            val = obj.combineBasis(obj.getDivDivBasis(pts,pts_ref,parts),parts,c);
        end

        function val = getD3Value(obj,pts,pts_ref,parts,c)
            %getD3Value - Evaluate third derivatives of the finite element function
            %   VAL = getD3Value(OBJ,PTS,PTS_REF,PARTS,C) evaluates the third
            %   derivatives of the finite element function with coefficients C
            %   in the elements PARTS at the reference point PTS_REF (1-by-d), or
            %   at the points PTS for elements that are not reference-element
            %   based, in the component order of getD3Basis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            curDoF4e = obj.dof4e(parts,:);
            u4e = reshape(c(curDoF4e),size(curDoF4e));
            d3 = FormBasis.elementRows(obj.getD3Basis(pts,pts_ref,parts),parts,obj);
            val = permute(sum(u4e.*d3,2),[1 3 2]);
        end

        function val = getLaplaceValue(obj,pts,pts_ref,parts,c)
            %getLaplaceValue - Evaluate the Laplacian of the finite element function
            %   VAL = getLaplaceValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the
            %   Laplacian of the finite element function with coefficients C in
            %   the elements PARTS at the reference point PTS_REF (1-by-d), or at
            %   the points PTS for elements that are not reference-element based,
            %   as the sum of the pure second derivatives of getD2Basis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            curDoF4e = obj.dof4e(parts,:);
            u4e = reshape(c(curDoF4e),size(curDoF4e));
            d2 = FormBasis.elementRows(obj.getD2Basis(pts,pts_ref,parts),parts,obj);
            val = sum(sum(u4e.*d2(:,:,1:size(obj.mesh.c4n,2)),3),2);
        end

        function val = getBiharmonicValue(obj,pts,pts_ref,parts,c)
            %getBiharmonicValue - Evaluate the biharmonic of the FE function
            %   VAL = getBiharmonicValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the
            %   biharmonic of the finite element function with coefficients C
            %   in the elements PARTS at the reference point PTS_REF (1-by-d), or
            %   at the points PTS for elements that are not reference-element
            %   based, as the sum of the fourth derivatives u_xxxx,u_yyyy(,u_zzzz)
            %   and twice the mixed u_xxyy(,u_xxzz,u_yyzz) of getD4Basis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            curDoF4e = obj.dof4e(parts,:);
            u4e = reshape(c(curDoF4e),size(curDoF4e));
            d4 = FormBasis.elementRows(obj.getD4Basis(pts,pts_ref,parts),parts,obj);
            dim = size(obj.mesh.c4n,2);
            if dim == 2
                mixed = 3;
            else
                mixed = 10:12;
            end
            val = sum(sum(u4e.*d4(:,:,1:dim),3),2) ...
                + 2*sum(sum(u4e.*d4(:,:,mixed),3),2);
        end

        function val = getTraceValue(obj,pts,pts_ref,parts,c)
            %getTraceValue - Evaluate the trace of the finite element function
            %   VAL = getTraceValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the
            %   trace of the matrix-valued finite element function with
            %   coefficients C in the elements PARTS at the reference point
            %   PTS_REF (1-by-d), or at the points PTS for elements that are not
            %   reference-element based, for elements that provide
            %   getTraceBasis.
            if ~isempty(pts_ref) && isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            curDoF4e = obj.dof4e(parts,:);
            u4e = reshape(c(curDoF4e),size(curDoF4e));
            traceBasis = FormBasis.elementRows(obj.getTraceBasis(pts,pts_ref,parts),parts,obj);
            val = sum(traceBasis.*u4e,2);
        end

        function x = interpolate(obj,func)
            %INTERPOLATE - Interpolate a function into the finite element space
            %   X = INTERPOLATE(OBJ,FUNC) evaluates FUNC at the points of
            %   the reference degrees of freedom c4refDof in all elements
            %   and assigns the values to the corresponding degrees of
            %   freedom. As in project, FUNC is a function of the physical
            %   points, @(x)f(x), or of the points, reference points and
            %   element indices, @(pts,pts_ref,parts)..., e.g. to
            %   interpolate another finite element function by its
            %   getValue (see FormEval.evalFunction). Elements without
            %   c4refDof raise an error; use project instead.
            x = zeros(obj.nrDoF,1);
            if isprop(obj,'c4refDof')
                parts = 1:obj.mesh.nrElems;
                for curRefDoF = 1 : size(obj.c4refDof,1)
                    pts_ref = obj.c4refDof(curRefDoF,:);
                    funcVal = FormEval.evalFunction(obj.mesh.point4e(pts_ref),pts_ref,parts,func,[]);
                    x(obj.dof4e(:,curRefDoF)) = funcVal;
                end
            else
                error('fem4m:FE:interpolationNotSupported','Interpolation not supported for this finite element.\n Use a projection instead.');
            end
        end

        function x = project(obj,func)
            %PROJECT - L2-project a function into the finite element space
            %   X = PROJECT(OBJ,FUNC) computes the coefficients of the local L2
            %   projections of the function handle FUNC by solving the local
            %   mass matrix systems on every element, scaled symmetrically
            %   by their diagonal, and averaging the values of shared
            %   degrees of freedom. In discontinuous spaces this is the L2
            %   projection; in the other spaces it reproduces the
            %   polynomials the space contains. FUNC is a function of the
            %   physical points, @(x)f(x), or of the points, reference
            %   points and element indices, @(pts,pts_ref,parts)..., e.g.
            %   to project another finite element function by its getValue.
            %
            %   Elements with elementwise degrees (allnrDoF) return all
            %   allnrDoF coefficients: the local problems use only the
            %   degrees of freedom used by the element, and the unused
            %   ones (nonDoFNodes) are zero.
            a = BilinearForm('u*v dx');
            [~,refMamaLoc] = a.assemble(obj);
            refMamaLoc = permute(refMamaLoc,[2 3 1]);
            b = LinearForm('f(x)*v dx');
            [~,b4e] = b.assemble(obj,'f',func);
            b4e = b4e';
            if isprop(obj,'allnrDoF')
                nrAllDoF = obj.allnrDoF;
            else
                nrAllDoF = obj.nrDoF;
            end
            isUsed = true(nrAllDoF,1);
            if isprop(obj,'nonDoFNodes')
                isUsed(obj.nonDoFNodes) = false;
            end
            x = zeros(nrAllDoF,1);
            for curElem = 1 : obj.mesh.nrElems
                dofs = obj.dof4e(curElem,:);
                used = isUsed(dofs);
                % scaled symmetrically by the diagonal: derivative dofs
                % (e.g. Argyris) scale the mass matrix by powers of h
                M = refMamaLoc(used,used,curElem);
                D = 1./sqrt(diag(M));
                x(dofs(used)) = x(dofs(used)) + D.*((M.*(D*D'))\(D.*b4e(used,curElem)));
            end
            count = accumarray(obj.dof4e(:),double(isUsed(obj.dof4e(:))),[nrAllDoF 1]);
            x(isUsed) = x(isUsed)./count(isUsed);
            x(~isUsed) = 0;
        end

        function P = prolongation(obj,feOld)
            %PROLONGATION - Prolongation matrix from a coarser space
            %   P = PROLONGATION(OBJ,FEOLD) returns the sparse matrix that maps
            %   the coefficients of a function of FEOLD, the same element
            %   on the mesh before the last refine call of OBJ.mesh, to
            %   the coefficients of the same function in OBJ. FEOLD may
            %   also be of a superclass of the class of OBJ with the same
            %   kind of values, e.g. HHJk for HHJB. P is exact for nested
            %   spaces and otherwise a local L2 projection (in reference
            %   coordinates for Piola-mapped spaces). The mesh provides
            %   only parent4e and its P1 prolongation; the local matrices
            %   come from the element itself (prolongationMatrix). Elements
            %   whose basis is not on the reference element (isRefElemBased
            %   false) or whose degrees of freedom are not numbered per
            %   element (DGPkFace) have no prolongation.
            P = prolongationMatrix(obj,feOld);
        end

        function x = prolongate(obj,xOld,feOld)
            %PROLONGATE - Prolongate coefficients of a coarser space
            %   X = PROLONGATE(OBJ,XOLD,FEOLD) maps the coefficients XOLD
            %   (one column per function) of FEOLD, the same element on
            %   the mesh before the last refine call of OBJ.mesh, to
            %   OBJ, i.e. X = prolongation(OBJ,FEOLD)*XOLD.
            x = obj.prolongation(feOld)*xOld;
        end

        function plot(obj,varargin)
            %PLOT - Plot the finite element function
            %   PLOT(OBJ,C) plots the finite element function with coefficients
            %   C by patches colored by its values at their vertices. On a 2d
            %   mesh it draws the graph of the function (on a surface mesh the
            %   colored surface), on a 3d mesh the Dirichlet, Neumann and
            %   symmetry boundary faces. Values with one component are drawn as
            %   they are; on a triangular mesh (also a surface mesh) values
            %   with three components [s11 s22 s12] or four components
            %   [s11 s12 s21 s22], on a 3d mesh values with six components
            %   [s11 s22 s33 s12 s13 s23] are colored by the von Mises
            %   stress, all other values by the Euclidean norm.
            %
            %   PLOT(OBJ,C,SUB) on a 2d mesh evaluates the function on the
            %   reference element refined uniformly SUB times (default 0).
            %
            %   PLOT(OBJ,C,SUB,...) on a 2d mesh passes the further arguments to
            %   patch.
            %
            %   PLOT(OBJ,C,X0,X1,Y0,Y1,Z0,Z1) on a 3d mesh draws the boundary of
            %   the elements whose four vertices lie in the box
            %   [X0,X1]x[Y0,Y1]x[Z0,Z1]; empty bounds are the bounds of the mesh,
            %   e.g. PLOT(OBJ,C,[],1/2,[],[],[],[]).
            %
            %   PLOT(OBJ,C,[],[],[],[],[],[],INDX) on a 3d mesh draws the
            %   boundary of the elements selected by INDX.
            if isa(obj.mesh,'Mesh')
                if numel(varargin) < 2
                    varargin{2} = [];
                end
                obj.plot2d(varargin{1},varargin{2},varargin{3:end});
            else
                if numel(varargin) < 8
                    varargin{8} = [];
                end
                obj.plot3d(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},varargin{7},varargin{8});
            end
        end

    end

    methods(Access=protected)

        function obj = dofManager(obj,orientEdges,orientFaces,byEntity)
            %dofManager - Enumerate the degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates, component after
            %   component, nrBasis4n nodal, nrBasis4ed edge, nrBasis4f face
            %   (3d only, zero if the element has no nrBasis4f), and
            %   nrBasis4e element degrees of freedom. The global numbering
            %   lists all nodes, then all edges, then all faces for each
            %   local index, followed by the element dofs element by
            %   element (so element-local dofs give block-diagonal element
            %   matrices), and dof4e lists the nodal, edge, face, and
            %   element dofs of each component in the same order. It sets
            %   dof4e and nrDoF and, if the element defines them, dof4ed and
            %   dof4f (in the global orientation of the edges and faces), and
            %   the dofs dof4Db, dof4Nb, and dof4Sb on the Dirichlet,
            %   Neumann, and symmetry boundary (the latter only if the mesh
            %   has n4Sb; its last column, the fixed direction, is
            %   ignored). The fixed and free degrees of freedom are left
            %   to the subclass.
            %
            %   OBJ = dofManager(OBJ,true) also reverses the order of the
            %   edge dofs in each element whose local edge orientation is
            %   opposite to the global one of n4ed (e.g. for point
            %   evaluations along the edge).
            %
            %   OBJ = dofManager(OBJ,ORIENTEDGES,true) also reorders the
            %   face dofs in each element by the permutation tables
            %   obj.fperm and obj.findices of the element: if the global
            %   node numbers of a local face are sorted by fperm(k,:),
            %   its i-th local face dof is the findices(k,i)-th dof of
            %   the face. The face dofs in dof4f are not reordered.
            %
            %   OBJ = dofManager(OBJ,ORIENTEDGES,ORIENTFACES,BYENTITY)
            %   with the logical vector BYENTITY = [NODES EDGES FACES]
            %   lists the local nodal, edge, or face dofs in dof4e entity
            %   by entity (all dofs of the first local entity, then of
            %   the second, ...) instead of index by index (the first
            %   dof of all local entities, then the second, ...). The
            %   global numbering, dof4ed, and dof4f do not change.
            if nargin < 2
                orientEdges = false;
            end
            if nargin < 3
                orientFaces = false;
            end
            if nargin < 4
                byEntity = false(1,3);
            end
            assert(isprop(obj,'nrBasis4n') && isprop(obj,'nrBasis4ed') && isprop(obj,'nrBasis4e'), ...
                'fem4m:FE:dofManager','Properties nrBasis4n, nrBasis4ed, or nrBasis4e not (all) defined.');
            assert(isequal(size(obj.nrBasis4n),size(obj.nrBasis4ed),size(obj.nrBasis4e)), ...
                'fem4m:FE:dofManager','The sizes of the dofs count per vertex, edge, and element do not match.');
            mesh = obj.mesh;
            hasFaces = isprop(mesh,'nrFaces');
            if isprop(obj,'nrBasis4f')
                nrBasis4f = obj.nrBasis4f;
            else
                nrBasis4f = zeros(size(obj.nrBasis4n));
            end
            assert(isequal(size(nrBasis4f),size(obj.nrBasis4n)), ...
                'fem4m:FE:dofManager','The sizes of the dofs count per vertex and face do not match.');
            assert(hasFaces || all(nrBasis4f == 0), ...
                'fem4m:FE:dofManager','Face degrees of freedom need a mesh with faces.');
            orientFaces = orientFaces && hasFaces;
            if orientFaces
                assert(isprop(obj,'fperm') && isprop(obj,'findices'), ...
                    'fem4m:FE:dofManager','Face orientation needs the properties fperm and findices.');
                assert(all(nrBasis4f == 0 | nrBasis4f == size(obj.findices,2)), ...
                    'fem4m:FE:dofManager','The columns of findices do not match nrBasis4f.');
                % node ids of the local faces (as in f4e) and their orientation
                n4f4e = [1 2 3; 1 2 4; 4 2 3; 1 3 4];
                perm4f = zeros(mesh.nrElems,4);
                for curF = 1 : 4
                    [~,sortIdx] = sort(mesh.n4e(:,n4f4e(curF,:)),2);
                    [~,perm4f(:,curF)] = ismember(sortIdx,obj.fperm,'rows');
                end
                assert(all(perm4f(:) > 0), ...
                    'fem4m:FE:dofManager','fperm does not contain all face node permutations.');
            end

            nrEd4e = size(mesh.ed4e,2);
            % first node of each local edge
            if nrEd4e == 3
                n4ed4e = [1 2 3];
            else
                n4ed4e = [1 2 3 1 2 3];
            end
            % optional arrays
            hasDof4ed = isprop(obj,'dof4ed');
            hasDof4f = isprop(obj,'dof4f') && hasFaces;
            hasDof4Db = isprop(obj,'dof4Db');
            hasDof4Nb = isprop(obj,'dof4Nb');
            hasDof4Sb = isprop(obj,'dof4Sb') && isprop(mesh,'n4Sb');
            if hasDof4ed; obj.dof4ed = zeros(mesh.nrEdges,0); end
            if hasDof4f; obj.dof4f = zeros(mesh.nrFaces,0); end
            if hasDof4Db; obj.dof4Db = zeros(0,1); end
            if hasDof4Nb; obj.dof4Nb = zeros(0,1); end
            if hasDof4Sb; obj.dof4Sb = zeros(0,1); end
            n4Db = unique(mesh.n4Db(:));
            n4Nb = unique(mesh.n4Nb(:));
            ed4Db = unique(mesh.ed4Db(:));
            ed4Nb = unique(mesh.ed4Nb(:));
            if hasFaces
                f4Db = unique(mesh.f4Db(:));
                f4Nb = unique(mesh.f4Nb(:));
            end
            if hasDof4Sb && ~isempty(mesh.n4Sb)
                n4Sb = unique(mesh.n4Sb(:,1:end-1));
                ed4Sb = unique(mesh.ed4Sb(:));
                f4Sb = zeros(0,1);
                if hasFaces; f4Sb = unique(mesh.f4Sb(:)); end
            else
                n4Sb = zeros(0,1);
                ed4Sb = zeros(0,1);
                f4Sb = zeros(0,1);
            end

            obj.dof4e = zeros(mesh.nrElems,0);
            nrDoF = 0;
            % component after component
            for component = 1 : numel(obj.nrBasis4n)
                % nodal dofs
                curDof4n = zeros(mesh.nrElems,0);
                for j = 1 : obj.nrBasis4n(component)
                    curDof4n = [curDof4n,nrDoF+mesh.n4e]; %#ok<AGROW>
                    if hasDof4ed; obj.dof4ed = [obj.dof4ed,nrDoF+mesh.n4ed]; end
                    if hasDof4f; obj.dof4f = [obj.dof4f,nrDoF+mesh.n4f]; end
                    if hasDof4Db; obj.dof4Db = [obj.dof4Db;nrDoF+n4Db]; end
                    if hasDof4Nb; obj.dof4Nb = [obj.dof4Nb;nrDoF+n4Nb]; end
                    if hasDof4Sb; obj.dof4Sb = [obj.dof4Sb;nrDoF+n4Sb]; end
                    nrDoF = nrDoF + mesh.nrNodes;
                end
                obj.dof4e = [obj.dof4e,entityMajor(curDof4n,byEntity(1),obj.nrBasis4n(component))];
                % edge dofs
                curDof4ed = zeros(mesh.nrElems,0);
                for j = 1 : obj.nrBasis4ed(component)
                    curDof4ed = [curDof4ed,nrDoF+mesh.ed4e]; %#ok<AGROW>
                    if hasDof4ed; obj.dof4ed = [obj.dof4ed,nrDoF+(1:mesh.nrEdges)']; end
                    if hasDof4f; obj.dof4f = [obj.dof4f,nrDoF+mesh.ed4f]; end
                    if hasDof4Db; obj.dof4Db = [obj.dof4Db;nrDoF+ed4Db]; end
                    if hasDof4Nb; obj.dof4Nb = [obj.dof4Nb;nrDoF+ed4Nb]; end
                    if hasDof4Sb; obj.dof4Sb = [obj.dof4Sb;nrDoF+ed4Sb]; end
                    nrDoF = nrDoF + mesh.nrEdges;
                end
                % fix edge orientation
                if orientEdges
                    for curEd = 1 : nrEd4e
                        permIdx = mesh.n4ed(mesh.ed4e(:,curEd),1) ~= mesh.n4e(:,n4ed4e(curEd));
                        curDof4ed(permIdx,curEd:nrEd4e:end) = ...
                            curDof4ed(permIdx,end-nrEd4e+curEd:-nrEd4e:curEd);
                    end
                end
                obj.dof4e = [obj.dof4e,entityMajor(curDof4ed,byEntity(2),obj.nrBasis4ed(component))];
                % face dofs
                curDof4f = zeros(mesh.nrElems,0);
                for j = 1 : nrBasis4f(component)
                    curDof4f = [curDof4f,nrDoF+mesh.f4e]; %#ok<AGROW>
                    if hasDof4f; obj.dof4f = [obj.dof4f,nrDoF+(1:mesh.nrFaces)']; end
                    if hasDof4Db; obj.dof4Db = [obj.dof4Db;nrDoF+f4Db]; end
                    if hasDof4Nb; obj.dof4Nb = [obj.dof4Nb;nrDoF+f4Nb]; end
                    if hasDof4Sb; obj.dof4Sb = [obj.dof4Sb;nrDoF+f4Sb]; end
                    nrDoF = nrDoF + mesh.nrFaces;
                end
                % fix face orientation
                if orientFaces && nrBasis4f(component) > 0
                    elems = (1:mesh.nrElems)'*ones(1,nrBasis4f(component));
                    for curF = 1 : 4
                        dof4f4e = curDof4f(:,curF:4:end);
                        curDof4f(:,curF:4:end) = dof4f4e(sub2ind(size(dof4f4e),elems,obj.findices(perm4f(:,curF),:)));
                    end
                end
                obj.dof4e = [obj.dof4e,entityMajor(curDof4f,byEntity(3),nrBasis4f(component))];
                % element dofs, element by element
                obj.dof4e = [obj.dof4e,nrDoF+reshape(1:obj.nrBasis4e(component)*mesh.nrElems,[],mesh.nrElems).'];
                nrDoF = nrDoF + obj.nrBasis4e(component)*mesh.nrElems;
            end
            obj.nrDoF = nrDoF;
        end

        function val = combineBasis(obj,basis,parts,c)
            %combineBasis - Combine evaluated basis functions with coefficients
            %   VAL = combineBasis(OBJ,BASIS,PARTS,C) sums the basis
            %   functions BASIS, evaluated in the elements PARTS, weighted
            %   with the coefficients C of the finite element function.
            curDoF4e = obj.dof4e(parts,:);
            coeff4e = reshape(c(curDoF4e),size(curDoF4e));
            basis = FormBasis.elementRows(basis,parts,obj);
            val = permute(sum(basis.*coeff4e,2),[1 3 2]);
        end

        function notDefined(obj,op)
            %notDefined - Raise the error for an undefined curl-type operator
            %   notDefined(OBJ,OP) raises the error fem4m:FE:<OP>NotDefined
            %   for the operator OP ('curl', 'rot', 'curlCurl', 'curlRot'
            %   or 'rotRot') that the element OBJ does not provide, with a
            %   hint at the operators defined for the dimension at hand.
            if isa(obj.mesh,'Mesh3d')
                first = 'in 3d, curl(u) applies to vector fields';
                hints = struct('curl',first, ...
                    'rot','rot(u) is the 2d rotation, use curl(u) in 3d', ...
                    'curlCurl','in 3d, curl(curl(u)) applies to vector fields', ...
                    'curlRot','curl(rot(u)) is a 2d operator, use curl(curl(u)) in 3d', ...
                    'rotRot','rot(rot(u)) is a 2d operator');
            else
                first = 'in 2d, curl(u) applies to scalar fields and rot(u) to vector and matrix fields';
                second = ['in 2d, curl(curl(u)) applies to scalar fields, curl(rot(u)) to vector ', ...
                    'fields and rot(rot(u)) to symmetric matrix fields'];
                hints = struct('curl',first,'rot',first, ...
                    'curlCurl',second,'curlRot',second,'rotRot',second);
            end
            names = struct('curl','curl(u)','rot','rot(u)','curlCurl','curl(curl(u))', ...
                'curlRot','curl(rot(u))','rotRot','rot(rot(u))');
            error(['fem4m:FE:',op,'NotDefined'],'%s: %s is not defined for this element; %s.', ...
                class(obj),names.(op),hints.(op));
        end

    end

    methods(Static,Access=protected)

        function val = mapDerivatives(c,d,order)
            %mapDerivatives - Map reference derivatives to the elements
            %   VAL = mapDerivatives(C,D,ORDER) maps the derivatives D of
            %   order ORDER on the reference element (one or nE rows, one
            %   column per function, the components of
            %   Monomials.multiIndices(n,ORDER) in the n reference
            %   variables) with the transformations C (nE x dim x n) of
            %   the elements, d_p = sum_q C(:,p,q) d_q. C may also hold
            %   the rows [c11,c21,c31,c12,...,c33] of
            %   Mesh3d.transformation. VAL holds the components of
            %   Monomials.multiIndices(dim,ORDER) on the elements.
            if ismatrix(c)
                c = permute(reshape(c,[],3,3),[1 3 2]);
            end
            [nE,dim,n] = size(c);
            alphaRef = Monomials.multiIndices(n,order);
            alphaPhys = Monomials.multiIndices(dim,order);
            % all ordered reference directions q, summed into their
            % components by B
            q = dec2base(0:n^order-1,n,order) - '0' + 1;
            [~,b] = ismember(reshape(sum(q == permute(1:n,[1 3 2]),2),[],n),alphaRef,'rows');
            B = double(b == 1:size(alphaRef,1));
            T = zeros(nE,size(alphaPhys,1),size(alphaRef,1));
            for a = 1 : size(alphaPhys,1)
                p = repelem(1:dim,alphaPhys(a,:));
                w = ones(nE,size(q,1));
                for s = 1 : order
                    w = w.*reshape(c(:,p(s),q(:,s)),nE,[]);
                end
                T(:,a,:) = reshape(w*B,nE,1,[]);
            end
            val = matMul(d,permute(T,[1 3 2]));
        end

        function val = curvedSurfaceD2(mesh,x,c,val,grad)
            %curvedSurfaceD2 - Correct the Hessians on a curved surface
            %   VAL = curvedSurfaceD2(MESH,X,C,VAL,GRAD) corrects the
            %   Hessians VAL [xx yy zz xy xz yz] that mapDerivatives maps
            %   from the reference triangle with the transformations
            %   C = MESH.transformation(X) of the curved surface MESH. It
            %   adds the derivative of C at the reference point X applied
            %   to the reference gradients GRAD (one or nE rows, one
            %   column per function) and projects the rows onto the
            %   tangent planes.
            val = reshape(val(:,:,[1 4 5 4 2 6 5 6 3]),size(val,1),size(val,2),3,3);
            gradT = permute(mesh.gradTransformation(x),[1 4 2 3]);
            for k = 1 : 3
                z = matMul(grad,gradT(:,:,:,k));
                z = matMul(z,permute(c,[1 3 2]));
                val(:,:,k,:) = val(:,:,k,:) + permute(z,[1 2 4 3]);
            end
            P = mesh.normalProjector(x);
            for i = 1 : size(val,2)
                val(:,i,:,:) = permute(matMul(P,squeeze(val(:,i,:,:))),[1 4 2 3]);
            end
            val = reshape(val,size(val,1),size(val,2),9);
            val = val(:,:,[1 5 9 2 3 6]);
        end

    end

    methods(Access=private)
        function plot2d(obj,c,sub,varargin)
            %plot2d - Plot the finite element function on a 2d mesh
            if nargin < 3 || isempty(sub)
                sub = 0;
            end
            meshRef = Mesh([0,0;1,0;0,1],[2 3 1],[1 2; 2 3; 3 1],[]);
            if sub > 0
                meshRef = meshRef.refine('uniform',sub);
            end
            color = zeros(obj.mesh.nrElems,meshRef.nrNodes);
            for j = 1: meshRef.nrNodes
                value = obj.getValue([],meshRef.c4n(j,:),[],c);
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
                    % norm
                    color(:,j) = sqrt(sum(value.^2,[2 3]));
                end
            end
            coordXref = reshape(meshRef.c4n(meshRef.n4e,1),[],3);
            coordYref = reshape(meshRef.c4n(meshRef.n4e,2),[],3);
            X = zeros(3,obj.mesh.nrElems);
            Y = zeros(3,obj.mesh.nrElems);
            if size(obj.mesh.c4n,2) == 3
                Z = zeros(3,obj.mesh.nrElems);
            end
            for j = 1 : meshRef.nrElems
                for i = 1 : 3
                    points = obj.mesh.point4e([coordXref(j,i)',coordYref(j,i)']);
                    X(i,:) = points(:,1);
                    Y(i,:) = points(:,2);
                    if size(obj.mesh.c4n,2) == 3
                        Z(i,:) = points(:,3);
                    end
                end
                C = color(:,meshRef.n4e(j,:))';
                if size(obj.mesh.c4n,2) == 3
                    patch(X,Y,Z,C,varargin{1:end});
                else
                    patch(X,Y,C,C,varargin{1:end});
                end
            end
        end

        function plot3d(obj,c,x0,x1,y0,y1,z0,z1,indx)
            %plot3d - Plot the finite element function on a 3d mesh
            if nargin <= 2 || (isempty(x0)&&isempty(x1)&&isempty(y0)&&isempty(y1)&&isempty(z0)&&isempty(z1)&&isempty(indx))
                n4fBd = [obj.mesh.n4Db;obj.mesh.n4Nb];
                if isprop(obj.mesh,'n4Sb') && ~isempty(obj.mesh.n4Sb)
                    n4fBd = [n4fBd;obj.mesh.n4Sb(:,1:3)];
                end
            else
                if nargin <= 8 || isempty(indx)
                    if nargin < 4 || isempty(x0); x0 = min(obj.mesh.c4n(:,1));end
                    if nargin < 5 || isempty(x1); x1 = max(obj.mesh.c4n(:,1));end
                    if nargin < 6 || isempty(y0); y0 = min(obj.mesh.c4n(:,2));end
                    if nargin < 7 || isempty(y1); y1 = max(obj.mesh.c4n(:,2));end
                    if nargin < 8 || isempty(z0); z0 = min(obj.mesh.c4n(:,3));end
                    if nargin < 9 || isempty(z1); z1 = max(obj.mesh.c4n(:,3));end
                    indx = x0 <= obj.mesh.c4n(obj.mesh.n4e(:,1),1) & obj.mesh.c4n(obj.mesh.n4e(:,1),1) <= x1 & ...
                        x0 <= obj.mesh.c4n(obj.mesh.n4e(:,2),1) & obj.mesh.c4n(obj.mesh.n4e(:,2),1) <= x1 & ...
                        x0 <= obj.mesh.c4n(obj.mesh.n4e(:,3),1) & obj.mesh.c4n(obj.mesh.n4e(:,3),1) <= x1 & ...
                        x0 <= obj.mesh.c4n(obj.mesh.n4e(:,4),1) & obj.mesh.c4n(obj.mesh.n4e(:,4),1) <= x1 & ...
                        y0 <= obj.mesh.c4n(obj.mesh.n4e(:,1),2) & obj.mesh.c4n(obj.mesh.n4e(:,1),2) <= y1 & ...
                        y0 <= obj.mesh.c4n(obj.mesh.n4e(:,2),2) & obj.mesh.c4n(obj.mesh.n4e(:,2),2) <= y1 & ...
                        y0 <= obj.mesh.c4n(obj.mesh.n4e(:,3),2) & obj.mesh.c4n(obj.mesh.n4e(:,3),2) <= y1 & ...
                        y0 <= obj.mesh.c4n(obj.mesh.n4e(:,4),2) & obj.mesh.c4n(obj.mesh.n4e(:,4),2) <= y1 & ...
                        z0 <= obj.mesh.c4n(obj.mesh.n4e(:,1),3) & obj.mesh.c4n(obj.mesh.n4e(:,1),3) <= z1 & ...
                        z0 <= obj.mesh.c4n(obj.mesh.n4e(:,2),3) & obj.mesh.c4n(obj.mesh.n4e(:,2),3) <= z1 & ...
                        z0 <= obj.mesh.c4n(obj.mesh.n4e(:,3),3) & obj.mesh.c4n(obj.mesh.n4e(:,3),3) <= z1 & ...
                        z0 <= obj.mesh.c4n(obj.mesh.n4e(:,4),3) & obj.mesh.c4n(obj.mesh.n4e(:,4),3) <= z1;
                end
                cutN4e = obj.mesh.n4e(indx,:);
                n4fBd = Mesh3d.getBoundary(cutN4e);
            end
            % sorted rows: the values below refer to the sorted face triples
            n4fBd = sort(n4fBd,2);
            [~,first,~,second] = intUniqueRows([n4fBd;obj.mesh.n4f]);
            f4Cut = second-size(n4fBd,1);
            f4Cut(second==first) = [];
            nodes = [0 0 ; 1 0 ; 0 1];
            color = zeros(numel(f4Cut),size(nodes,1));
            for j = 1 : size(nodes,1)
                value = getFaceValues4RefElemFunktion(nodes(j,:),@(z)obj.getValue([],z,[],c),obj.mesh);
                value = value(f4Cut,:,1);
                if size(value,2) == 1
                    % value
                    color(:,j) = value;
                elseif size(value,2) == 6
                    % von Mises stress
                    color(:,j) = sqrt(1/2*((value(:,1)-value(:,2)).^2+(value(:,2)-value(:,3)).^2+(value(:,3)-value(:,1)).^2) ...
                        +3*(value(:,4).^2+value(:,5).^2+value(:,6).^2));
                else
                    % norm
                    color(:,j) = sqrt(sum(value.^2,[2 3]));
                end
            end
            X = obj.mesh.c4n(:,1);
            Y = obj.mesh.c4n(:,2);
            Z = obj.mesh.c4n(:,3);
            patch(X(n4fBd)',Y(n4fBd)',Z(n4fBd)',color');
        end

    end

end

function dofs = entityMajor(dofs,byEntity,k)
%entityMajor - Reorder a block of local dofs entity by entity
%   DOFS = entityMajor(DOFS,BYENTITY,K) reorders the columns of DOFS,
%   given index by index (the first dof of all local entities, then the
%   second, ...) with K dofs per entity, entity by entity (all K dofs of
%   the first entity, then of the second, ...) if BYENTITY is true.
if ~byEntity || k <= 1
    return
end
nrEntities = size(dofs,2)/k;
dofs = dofs(:,reshape(reshape(1:nrEntities*k,nrEntities,k)',1,[]));
end
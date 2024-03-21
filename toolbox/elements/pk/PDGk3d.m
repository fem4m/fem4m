classdef PDGk3d < FE
    %PDGk3d - Discontinuous scaled monomial element of degree k on tetrahedra
    %   OBJ = PDGk3d(DEGREE,MESH) creates the discontinuous
    %   element of polynomial degree DEGREE >= 0 on the tetrahedral
    %   mesh MESH (Mesh3d). Its basis functions on an element are the
    %   monomials of total degree at most DEGREE in the element-scaled
    %   coordinates (x-mid4e)/diam4e, with the centroid mid4e and the
    %   diameter diam4e of the element, in the order of
    %   Monomials.exponents.
    %
    %   The basis is evaluated at the physical points. Without
    %   physical points it is evaluated at the images of the
    %   reference point in the elements.
    %
    %   PDGk3d functions:
    %       getBasis     - Evaluate the basis functions
    %       getGradBasis - Evaluate the gradients of the basis
    %       getD2Basis   - Evaluate the Hessians of the basis
    %       getD3Basis   - Evaluate third derivatives of the basis
    %       getD4Basis   - Evaluate fourth derivatives of the basis
    %
    %   PDGk3d properties:
    %       mesh           - Mesh on which the element is defined
    %       nrDoF          - Total number of degrees of freedom
    %       dof4e          - Local-to-global dofs per element
    %       freeNodes      - Indices of the free degrees of freedom
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       degree         - Polynomial degree of the element
    %       nrBasisFunc    - Number of basis functions per element
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of edge basis functions
    %       nrBasis4f      - Number of face basis functions
    %       nrBasis4e      - Number of interior basis functions
    %       isRefElemBased - Basis on the physical elements
    %
    %   See also PDGk, DGDk3d, Pk3d, Monomials, Mesh3d

    % Copyright (c) 2026 by the fem4m authors

    properties (Access = public)
        mesh               % mesh on which the element is defined
        nrDoF              % total number of degrees of freedom
        dof4e              % local-to-global dofs per element
        freeNodes          % indices of the free degrees of freedom
        fixedNodes         % indices of the fixed degrees of freedom
        degree             % polynomial degree of the element
        nrBasisFunc        % number of basis functions per element
        nrBasis4n          % number of nodal basis functions
        nrBasis4ed         % number of edge basis functions
        nrBasis4f          % number of face basis functions
        nrBasis4e          % number of interior basis functions
        isRefElemBased = false % basis on the physical elements
    end

    properties (Access = private)
        monomBasis             % handle of the monomials, one row per point
        gradBasisCoefficients  % gradient coefficients of the monomials
        d2BasisCoefficients    % Hessian coefficients of the monomials
        d3BasisCoefficients    % third derivative coefficients
        d4BasisCoefficients    % fourth derivative coefficients
    end

    methods
        function obj = PDGk3d(degree,mesh)
            %PDGk3d - Create a PDGk3d finite element object
            if degree < 0
                error('fem4m:PDGk3d:invalidDegree','PDGk3d degree must be at least 0');
            end
            obj.degree = degree;
            obj.mesh = mesh;
            obj = obj.getMonomBasis;
            obj = obj.dofCounts();
            obj = obj.dofManager;
        end

        function val = getBasis(obj,pts,pts_ref,parts)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   scaled monomial basis at the physical points PTS
            %   (N-by-3), row i in the element PARTS(i), one row per
            %   point. If PTS is empty, the points are the images of
            %   the reference point PTS_REF (1-by-3) in the elements
            %   PARTS. An empty PARTS means all elements.
            val = obj.getScaledMonomBasis(pts,pts_ref,parts);
        end

        function val = getGradBasis(obj,pts,pts_ref,parts)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,PTS_REF,PARTS) evaluates
            %   the gradients [u_x,u_y,u_z] of the scaled monomial
            %   basis at the points of getBasis. VAL has one row per
            %   point, one column per basis function and the components
            %   along the third dimension.
            val = obj.getScaledMonomDerivBasis(pts,pts_ref,parts,obj.gradBasisCoefficients,1);
        end

        function val = getD2Basis(obj,pts,pts_ref,parts)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   second derivatives [u_xx,u_yy,u_zz,u_xy,u_xz,u_yz] of
            %   the scaled monomial basis at the points of getBasis.
            %   VAL is laid out as in getGradBasis.
            val = obj.getScaledMonomDerivBasis(pts,pts_ref,parts,obj.d2BasisCoefficients,2);
        end

        function val = getD3Basis(obj,pts,pts_ref,parts)
            %getD3Basis - Evaluate third derivatives of the basis
            %   VAL = getD3Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   ten third derivatives [xxx,yyy,zzz,xxy,xxz,xyy,yyz,
            %   xzz,yzz,xyz] of the scaled monomial basis at the
            %   points of getBasis.
            %   VAL is laid out as in getGradBasis.
            val = obj.getScaledMonomDerivBasis(pts,pts_ref,parts,obj.d3BasisCoefficients,3);
        end

        function val = getD4Basis(obj,pts,pts_ref,parts)
            %getD4Basis - Evaluate fourth derivatives of the basis
            %   VAL = getD4Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   15 fourth derivatives [xxxx,yyyy,zzzz,xxxy,xxxz,xyyy,
            %   yyyz,xzzz,yzzz,xxyy,xxzz,yyzz,xxyz,xyyz,xyzz] of the
            %   scaled monomial basis at the points of getBasis.
            %   VAL is laid out as in getGradBasis.
            val = obj.getScaledMonomDerivBasis(pts,pts_ref,parts,obj.d4BasisCoefficients,4);
        end

    end

    methods (Access = private)

        function [z,h] = scaledCoordinates(obj,pts,pts_ref,parts)
            %scaledCoordinates - Element-scaled coordinates of the points
            %   [Z,H] = scaledCoordinates(OBJ,PTS,PTS_REF,PARTS) returns
            %   the coordinates Z = (PTS-mid4e)/diam4e of the points of
            %   getBasis in the elements PARTS and the diameters H of
            %   these elements.
            if isempty(parts)
                parts = 1:obj.mesh.nrElems;
            end
            if isempty(pts)
                pts = obj.mesh.point4e(pts_ref);
                pts = pts(parts,:);
            end
            mid4e = obj.mesh.mid4e;
            diam4e = obj.mesh.diam4e;
            h = diam4e(parts);
            z = (pts-mid4e(parts,:))./h(:);
        end

        function val = getScaledMonomBasis(obj,pts,pts_ref,parts)
            %getScaledMonomBasis - Evaluate the scaled monomial basis
            z = obj.scaledCoordinates(pts,pts_ref,parts);
            val = obj.monomBasis(z(:,1),z(:,2),z(:,3));
        end

        function val = getScaledMonomDerivBasis(obj,pts,pts_ref,parts,coeff,order)
            %getScaledMonomDerivBasis - Evaluate scaled monomial derivatives
            %   VAL = getScaledMonomDerivBasis(OBJ,PTS,PTS_REF,PARTS,COEFF,ORDER)
            %   evaluates the derivatives of order ORDER of the scaled
            %   monomials with the stacked coefficients COEFF at the
            %   points of getBasis, one component per block of COEFF,
            %   times diam4e^-ORDER by the chain rule.
            [z,h] = obj.scaledCoordinates(pts,pts_ref,parts);
            m = obj.monomBasis(z(:,1),z(:,2),z(:,3));
            val = reshape(m*coeff.',size(m,1),size(m,2),[])./h(:).^order;
        end

        function obj = getMonomBasis(obj)
            %getMonomBasis - Build the monomial basis functions
            %   OBJ = getMonomBasis(OBJ) sets nrBasisFunc, the handle
            %   monomBasis of the monomials x^i*y^j*z^k of total degree
            %   up to obj.degree in the lexicographic order of
            %   Monomials.exponents, and the matrices
            %   gradBasisCoefficients, d2BasisCoefficients,
            %   d3BasisCoefficients, and d4BasisCoefficients of their
            %   first to fourth derivatives.
            E = Monomials.exponents(3,obj.degree);
            obj.nrBasisFunc = size(E,1);
            obj.monomBasis = Monomials.basis(E);
            obj.gradBasisCoefficients = Monomials.derivatives(E,Monomials.multiIndices(3,1));
            obj.d2BasisCoefficients = Monomials.derivatives(E,Monomials.multiIndices(3,2));
            obj.d3BasisCoefficients = Monomials.derivatives(E,Monomials.multiIndices(3,3));
            obj.d4BasisCoefficients = Monomials.derivatives(E,Monomials.multiIndices(3,4));
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, nrBasis4f, and nrBasis4e of nodal, edge,
            %   face, and interior degrees of freedom. The constructor
            %   calls it right before dofManager, so that subclasses may
            %   override it.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 0;
            obj.nrBasis4f = 0;
            obj.nrBasis4e = obj.nrBasisFunc;
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) numbers the element-local
            %   degrees of freedom with FE.dofManager and keeps all of
            %   them free.
            obj = dofManager@FE(obj);
            obj.fixedNodes = [];
            obj.freeNodes = 1:obj.nrDoF;
        end

    end

end
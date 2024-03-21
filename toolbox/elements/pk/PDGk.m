classdef PDGk < FE
    %PDGk - Discontinuous element of scaled monomials of degree k
    %   OBJ = PDGk(DEGREE,MESH) creates the discontinuous element of
    %   polynomial degree DEGREE >= 0 on the polygonal mesh MESH
    %   (MeshPoly). Its basis functions on an element are the
    %   monomials of total degree at most DEGREE in the element-scaled
    %   coordinates (x-mid4e)/diam4e, with the centroid mid4e and the
    %   diameter diam4e of the element, in the order of
    %   Monomials.exponents. The basis lives on the physical elements
    %   and is evaluated at physical points.
    %
    %   PDGk functions:
    %       getBasis         - Evaluate the basis functions
    %       getGradBasis     - Evaluate the gradients of the basis
    %       getD2Basis       - Evaluate the Hessians of the basis
    %       getD3Basis       - Evaluate third derivatives of the basis
    %       getD4Basis       - Evaluate fourth derivatives of the basis
    %       plot             - Plot the finite element function
    %       projectSubmeshDG - Project onto a DG subtriangulation
    %
    %   PDGk properties:
    %       mesh           - Mesh on which the element is defined
    %       nrDoF          - Total number of degrees of freedom
    %       dof4e          - Local-to-global dofs per element
    %       freeNodes      - Indices of the free degrees of freedom
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       degree         - Polynomial degree of the element
    %       nrBasisFunc    - Number of basis functions per element
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of edge basis functions
    %       nrBasis4e      - Number of interior basis functions
    %       isRefElemBased - Basis on the physical elements
    %
    %   See also Pk, DGDk, PDGk3d, PDGkVector, PDGkMatrix, MeshPoly

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
        function obj = PDGk(degree,mesh)
            %PDGk - Create a PDGk finite element object
            if degree < 0
                error('fem4m:PDGk:invalidDegree','PDGk degree must be at least 0');
            end
            obj.degree = degree;
            obj.mesh = mesh;
            obj = obj.getMonomBasis;
            obj = obj.dofCounts();
            obj = obj.dofManager;
        end

        function val = getBasis(obj,pts,~,parts)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   scaled monomial basis at the physical points PTS
            %   (N-by-2), row i in the element PARTS(i), one row per
            %   point. PTS_REF is ignored.
            val = obj.getScaledMonomBasis(pts,parts);
        end

        function val = getGradBasis(obj,pts,~,parts)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   gradients [u_x,u_y] of the scaled monomial basis at the
            %   physical points PTS in the elements PARTS as in
            %   getBasis. PTS_REF is ignored. VAL has one row per point,
            %   one column per basis function and the components along
            %   the third dimension.
            val = obj.getScaledMonomDerivBasis(pts,parts,obj.gradBasisCoefficients,1);
        end

        function val = getD2Basis(obj,pts,~,parts)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   second derivatives [u_xx,u_yy,u_xy] of the scaled
            %   monomial basis at the physical points PTS in the
            %   elements PARTS as in getBasis. PTS_REF is ignored.
            %   VAL is laid out as in getGradBasis.
            val = obj.getScaledMonomDerivBasis(pts,parts,obj.d2BasisCoefficients,2);
        end

        function val = getD3Basis(obj,pts,~,parts)
            %getD3Basis - Evaluate third derivatives of the basis
            %   VAL = getD3Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   third derivatives [u_xxx,u_yyy,u_xxy,u_xyy] of the
            %   scaled monomial basis at the physical points PTS in the
            %   elements PARTS as in getBasis. PTS_REF is ignored.
            %   VAL is laid out as in getGradBasis.
            val = obj.getScaledMonomDerivBasis(pts,parts,obj.d3BasisCoefficients,3);
        end

        function val = getD4Basis(obj,pts,~,parts)
            %getD4Basis - Evaluate fourth derivatives of the basis
            %   VAL = getD4Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   fourth derivatives [u_xxxx,u_yyyy,u_xxyy,u_xxxy,u_xyyy]
            %   of the scaled monomial basis at the physical points PTS
            %   in the elements PARTS as in getBasis. PTS_REF is ignored.
            %   VAL is laid out as in getGradBasis.
            val = obj.getScaledMonomDerivBasis(pts,parts,obj.d4BasisCoefficients,4);
        end

        function plot(obj,c)
            %PLOT - Plot the finite element function
            %   PLOT(OBJ,C) plots the piecewise polynomial function
            %   with coefficients C over the polygonal mesh: every
            %   polygon is drawn as a patch whose vertices are lifted to
            %   the values of the function there and colored by them,
            %   interpolated in between. The axes are held (hold on).
            X = cellfun(@(curElem)(obj.mesh.c4n(curElem,1)),obj.mesh.n4e,'UniformOutput',false)';
            Y = cellfun(@(curElem)(obj.mesh.c4n(curElem,2)),obj.mesh.n4e,'UniformOutput',false)';
            for j = unique(obj.mesh.nrN4e)'
                parts = find(obj.mesh.nrN4e==j);
                curX = cell2mat(X(1,parts));
                curY = cell2mat(Y(1,parts));
                curZ = zeros(size(curY));
                for k = 1 : j
                    curZ(k,:) = obj.getValue([curX(k,:)',curY(k,:)'],[],parts,c);
                end
                patch(curX,curY,curZ,curZ);
                hold on;
            end
        end

        function [c,dgfe] = projectSubmeshDG(obj,c,dgfe)
            %projectSubmeshDG - Project onto a DG subtriangulation
            %   C = projectSubmeshDG(OBJ,C) projects the piecewise
            %   polynomial with coefficients C onto the DGDk element of
            %   the same degree on the centroid subtriangulation
            %   subtriangulateMid of the mesh and returns the
            %   coefficients there. The projection solves the discrete
            %   problem of the form 'grad(u)*grad(v)dx + h^-2*u*v dx'
            %   with the right-hand side given by the function and its
            %   gradient.
            %
            %   [C,DGFE] = projectSubmeshDG(OBJ,C) also returns
            %   the element DGFE used for the projection.
            %
            %   [...] = projectSubmeshDG(OBJ,C,DGFE) projects
            %   onto the given element DGFE, whose mesh is a
            %   subtriangulation of the mesh with its parent4e set,
            %   e.g. from subtriangulateMid. An empty DGFE builds the
            %   default one.
            if nargin < 3 || isempty(dgfe)
                submesh = obj.mesh.subtriangulateMid;
                dgfe = DGDk(obj.degree,submesh);
            end
            parent4e = dgfe.mesh.parent4e;
            % project onto subtriangulation
            m = BilinearForm('grad(u)*grad(v)dx + h^-2*u*v dx');
            M = m.assemble(dgfe);
            b = LinearForm('h^-2*f(x)*v dx + g(x)*grad(v) dx');
            rhs = b.assemble(dgfe,'f',@(pts,pts_ref,ind)(obj.getValue(pts,[],parent4e(ind),c)),...
                'g',@(pts,pts_ref,ind)(obj.getGradValue(pts,[],parent4e(ind),c)));
            c = M\rhs;
        end

    end

    methods (Access = private)

        function val = getScaledMonomBasis(obj,pts,parts)
            %getScaledMonomBasis - Evaluate the scaled monomial basis
            mid4e = obj.mesh.mid4e;
            diam4e = obj.mesh.diam4e;
            z = (pts-mid4e(parts,:))./diam4e(parts);
            val = obj.monomBasis(z(:,1),z(:,2));
        end

        function val = getScaledMonomDerivBasis(obj,pts,parts,coeff,order)
            %getScaledMonomDerivBasis - Evaluate scaled monomial derivatives
            %   VAL = getScaledMonomDerivBasis(OBJ,PTS,PARTS,COEFF,ORDER)
            %   evaluates the derivatives of order ORDER of the scaled
            %   monomials with the stacked coefficients COEFF at the
            %   physical points PTS in the elements PARTS, one component
            %   per block of COEFF, times diam4e^-ORDER by the chain rule.
            mid4e = obj.mesh.mid4e;
            diam4e = obj.mesh.diam4e;
            z = (pts-mid4e(parts,:))./diam4e(parts);
            m = obj.monomBasis(z(:,1),z(:,2));
            val = reshape(m*coeff.',size(m,1),size(m,2),[])./diam4e(parts).^order;
        end

        function obj = getMonomBasis(obj)
            %getMonomBasis - Build the monomial basis functions
            %   OBJ = getMonomBasis(OBJ) sets nrBasisFunc, the handle
            %   monomBasis of the monomials x^i*y^j of total degree up
            %   to obj.degree in the lexicographic order of
            %   Monomials.exponents, and the matrices gradBasisCoefficients,
            %   d2BasisCoefficients, d3BasisCoefficients, and
            %   d4BasisCoefficients of their first to fourth derivatives.
            E = Monomials.exponents(2,obj.degree);
            obj.nrBasisFunc = size(E,1);
            obj.monomBasis = Monomials.basis(E);
            obj.gradBasisCoefficients = Monomials.derivatives(E,Monomials.multiIndices(2,1));
            obj.d2BasisCoefficients = Monomials.derivatives(E,Monomials.multiIndices(2,2));
            obj.d3BasisCoefficients = Monomials.derivatives(E,Monomials.multiIndices(2,3));
            obj.d4BasisCoefficients = Monomials.derivatives(E,Monomials.multiIndices(2,4));
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of nodal, edge, and interior
            %   degrees of freedom. The constructor calls it right
            %   before dofManager, so that subclasses may override it.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 0;
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
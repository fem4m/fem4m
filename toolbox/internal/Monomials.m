classdef Monomials
    %MONOMIALS - Monomials and the matrices of their derivatives
    %   Monomials provides static methods for the monomials given by an
    %   exponent table E with one row [i j] or [i j k] per monomial
    %   x^i*y^j or x^i*y^j*z^k. The partial derivatives of the column m
    %   of the monomials are again combinations of the monomials,
    %   d^alpha m = D*m, with small integer matrices D, so that the
    %   derivatives of a basis C*m of any order follow from the one
    %   field of monomial values as (C*D)*m. The elements use the
    %   lexicographic exponent tables returned by exponents, and
    %   massMatrix integrates the products of the monomials exactly.
    %
    %   MONOMIALS functions:
    %       exponents         - Exponents of the monomials up to a degree
    %       basis             - Function handle of the monomials
    %       multiIndices      - Multi-indices of the derivatives of an order
    %       derivatives       - Derivative matrices of the monomials
    %       basisCoefficients - Derivative coefficients of a basis
    %       massMatrix        - Integrals of the monomial products
    %
    %   See also Pk, Pk3d, PDGk, BarycentricJacobi, JacobiPolynomials

    % Copyright (c) 2026 by the fem4m authors

    methods (Static)

        function E = exponents(dim,degree)
            %EXPONENTS - Exponents of the monomials up to a degree
            %   E = EXPONENTS(DIM,DEGREE) returns one row [i j] (DIM = 2)
            %   or [i j k] (DIM = 3) per monomial of total degree at
            %   most DEGREE, in lexicographic order (first exponent
            %   outer, last inner).
            g = cell(1,dim);
            [g{:}] = ndgrid(0:degree);
            E = reshape(cat(dim+1,g{:}),[],dim);
            E = sortrows(E(sum(E,2) <= degree,:));
        end

        function f = basis(E)
            %BASIS - Function handle of the monomials
            %   F = BASIS(E) returns the handle F(X,Y), or F(X,Y,Z) for
            %   three columns of E, of the monomials with the exponents
            %   E. The coordinates are vectors of equal length or
            %   scalars, and F returns one row per point and one column
            %   per monomial.
            ex = E(:,1)';
            ey = E(:,2)';
            if size(E,2) == 2
                f = @(x,y)x(:).^ex.*y(:).^ey;
            else
                ez = E(:,3)';
                f = @(x,y,z)x(:).^ex.*y(:).^ey.*z(:).^ez;
            end
        end

        function alpha = multiIndices(dim,order)
            %multiIndices - Multi-indices of the derivatives of an order
            %   ALPHA = multiIndices(DIM,ORDER) returns the multi-indices
            %   of the partial derivatives of order ORDER = 1,...,4 in
            %   DIM = 2 or 3 variables, one row per component, in the
            %   component order of the elements: [x y], [xx yy xy],
            %   [xxx yyy xxy xyy], and [xxxx yyyy xxyy xxxy xyyy] in 2d,
            %   [x y z], [xx yy zz xy xz yz], [xxx yyy zzz xxy xxz xyy
            %   yyz xzz yzz xyz], and [xxxx yyyy zzzz xxxy xxxz xyyy yyyz
            %   xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz] in 3d.
            if dim == 2
                tables = {[1 0;0 1],[2 0;0 2;1 1],[3 0;0 3;2 1;1 2], ...
                    [4 0;0 4;2 2;3 1;1 3]};
            else
                tables = {eye(3),[2 0 0;0 2 0;0 0 2;1 1 0;1 0 1;0 1 1], ...
                    [3 0 0;0 3 0;0 0 3;2 1 0;2 0 1;1 2 0;0 2 1;1 0 2;0 1 2;1 1 1], ...
                    [4 0 0;0 4 0;0 0 4;3 1 0;3 0 1;1 3 0;0 3 1;1 0 3;0 1 3; ...
                    2 2 0;2 0 2;0 2 2;2 1 1;1 2 1;1 1 2]};
            end
            alpha = tables{order};
        end

        function D = derivatives(E,alpha)
            %DERIVATIVES - Derivative matrices of the monomials
            %   D = DERIVATIVES(E,ALPHA) returns D = [D_1;...;D_K] with
            %   d^ALPHA(c,:) m = D_c*m for the column m of the monomials
            %   with the exponents E and the K = size(ALPHA,1)
            %   multi-indices ALPHA, one per row. D_c is the product of
            %   the powers D_r^ALPHA(c,r) of the first derivative
            %   matrices D_r, whose row a holds the exponent E(a,r) in
            %   the column of the monomial with the exponents
            %   E(a,:)-e_r. Each D_c is n-by-n for the n = size(E,1)
            %   monomials. Derivatives of monomials outside E are
            %   dropped, so E must contain E(a,:)-e_r for every row a
            %   with E(a,r) > 0, as the tables of exponents do.
            [n,dim] = size(E);
            Dr = cell(1,dim);
            for r = 1 : dim
                F = E;
                F(:,r) = F(:,r)-1;
                [found,col] = ismember(F,E,'rows');
                row = find(found);
                Dr{r} = zeros(n,n);
                Dr{r}(sub2ind([n n],row,col(row))) = E(row,r);
            end
            D = zeros(n*size(alpha,1),n);
            for c = 1 : size(alpha,1)
                Dc = eye(n);
                for r = 1 : dim
                    Dc = Dc*Dr{r}^alpha(c,r);
                end
                D((c-1)*n+(1:n),:) = Dc;
            end
        end

        function G = basisCoefficients(C,E,order)
            %basisCoefficients - Derivative coefficients of a basis
            %   G = basisCoefficients(C,E,ORDER) returns the stacked
            %   coefficients G = [C*D_1;...;C*D_K] of the derivatives of
            %   order ORDER of the basis C*m, m the column of the
            %   monomials with the exponents E, D_c from derivatives and
            %   the K components from multiIndices. C has one row per
            %   basis function and one column per monomial. Block c of
            %   G*m holds component c of the derivatives of all basis
            %   functions.
            D = Monomials.derivatives(E,Monomials.multiIndices(size(E,2),order));
            n = size(E,1);
            G = reshape(C*reshape(D,n,[]),[],n);
        end

        function M = massMatrix(E)
            %massMatrix - Integrals of the monomial products
            %   M = massMatrix(E) returns the integrals
            %   M(a,b) = alpha!/(|alpha|+d)!, alpha = E(a,:)+E(b,:), of
            %   the products of the monomials with the exponents E over
            %   the reference simplex conv{0,e_1,...,e_d} in
            %   d = size(E,2) dimensions.
            S = permute(E,[1 3 2]) + permute(E,[3 1 2]);
            M = prod(factorial(S),3)./factorial(sum(S,3)+size(E,2));
        end

    end

end
classdef BarycentricJacobi
    %BarycentricJacobi - Derivatives of barycentric Jacobi factors
    %   BarycentricJacobi provides static methods for the factors
    %   L_j(lambda_b,sum(lambda_nodes)), j = 0,...,N, of the
    %   hierarchical bases of Beuchler, Pillwein, Schoeberl, and
    %   Zaglmayr on the reference triangle or tetrahedron, where
    %   L_j(x,s) = s^j L_j^ALPHA(2x/s-1) are the scaled integrated
    %   Jacobi polynomials of JacobiPolynomials.scaledIntegratedPoly
    %   and lambda the barycentric coordinates. The factors and their
    %   products are polynomials, also on the edges and at the
    %   vertices. The derivatives of all orders up to M <= 4 are
    %   stored as columns in the order of the rows of
    %   multiIndices(DIM,M): the value, the gradient, and the
    %   components of Monomials.multiIndices of the higher orders.
    %
    %   BarycentricJacobi functions:
    %       coordinates  - Barycentric coordinates and their gradients
    %       chainRule    - Chain rule for the arguments of a factor
    %       factor       - Derivatives of a scaled integrated Jacobi factor
    %       leibnizTerms - Terms of the Leibniz rule of a product
    %       product      - Derivatives of products of two factors
    %       whitney      - Whitney edge function and its curl
    %       whitneyFace  - Whitney face function and its divergence
    %       timesField   - Products of functions with a vector field
    %       cross        - Cross products of rows
    %       multiIndices - Multi-indices of all derivatives up to an order
    %       columns      - Columns of the derivatives of one order
    %
    %   See also JacobiPolynomials, Monomials, HPk3d, Nd1k3d, RTk3d

    % Copyright (c) 2026 by the fem4m authors

    methods (Static)

        function [lambda,gradLambda] = coordinates(x)
            %COORDINATES - Barycentric coordinates and their gradients
            %   LAMBDA = COORDINATES(X) returns the barycentric
            %   coordinates LAMBDA = [1-sum(X); X(:)], a (d+1)-by-1
            %   column, of the point X (1-by-d) of the reference triangle
            %   (d = 2) or tetrahedron (d = 3).
            %
            %   [LAMBDA,GRADLAMBDA] = COORDINATES(X) also returns the
            %   constant gradients GRADLAMBDA of the coordinates, a
            %   (d+1)-by-d matrix with one row per coordinate.
            dim = numel(x);
            lambda = [1-sum(x); x(:)];
            gradLambda = [-ones(1,dim); eye(dim)];
        end

        function chain = chainRule(b,nodes,m,dim)
            %chainRule - Chain rule for the arguments of a factor
            %   CHAIN = chainRule(B,NODES,M,DIM) returns the chain rule
            %   of the derivatives up to order M <= 4 of a function
            %   f(X,S) of X = lambda_B and S = sum(lambda_NODES) in
            %   DIM = 2 or 3 variables: the derivatives of
            %   multiIndices(DIM,M) are CHAIN.C times the derivatives
            %   d^kx/dX^kx d^ks/dS^ks f with kx = CHAIN.kx(c) and
            %   ks = CHAIN.ks(c) for the column c of CHAIN.C, ordered by
            %   their order and, within an order, by decreasing kx.
            %   CHAIN also holds B and NODES in the fields b and nodes.
            [~,gradLambda] = BarycentricJacobi.coordinates(zeros(1,dim));
            gX = gradLambda(b,:);
            gS = sum(gradLambda(nodes,:),1);
            D = BarycentricJacobi.multiIndices(dim,m);
            key = (m+1).^(0:dim-1)';
            row4D = zeros((m+1)^dim,1);
            row4D(D*key+1) = 1:size(D,1);
            C = zeros(size(D,1),(m+1)*(m+2)/2);
            C(1,1) = 1;
            % d^D = (gX(r) d/dX + gS(r) d/dS) d^(D-e_r)
            for l = 2 : size(D,1)
                r = find(D(l,:),1);
                prev = D(l,:);
                prev(r) = prev(r) - 1;
                k = sum(D(l,:));
                C(l,k*(k+1)/2+(1:k+1)) = conv(C(row4D(prev*key+1),(k-1)*k/2+(1:k)),[gX(r) gS(r)]);
            end
            k = repelem(0:m,1:m+1);
            ks = (0:numel(k)-1) - k.*(k+1)/2;
            chain = struct('b',b,'nodes',nodes,'C',C,'kx',k-ks,'ks',ks);
        end

        function F = factor(lambda,chain,n,alpha)
            %FACTOR - Derivatives of a scaled integrated Jacobi factor
            %   F = FACTOR(LAMBDA,CHAIN,N,ALPHA) evaluates the
            %   derivatives of the factors L_j^ALPHA(lambda_b,
            %   sum(lambda_nodes)), j = 0,...,N, at the barycentric
            %   coordinates LAMBDA of one point with b, nodes, and the
            %   chain rule of CHAIN from chainRule. F has one row per
            %   degree j and one column per derivative of multiIndices;
            %   for N < 0 it has no rows.
            F = zeros(max(n+1,0),size(chain.C,1));
            if n < 0
                return;
            end
            X = lambda(chain.b);
            S = sum(lambda(chain.nodes));
            for c = find(any(chain.C,1))
                F = F + JacobiPolynomials.scaledIntegratedPoly(X,S,n,alpha,chain.kx(c),chain.ks(c))*chain.C(:,c)';
            end
        end

        function terms = leibnizTerms(D,m)
            %leibnizTerms - Terms of the Leibniz rule of a product
            %   TERMS = leibnizTerms(D,M) returns the terms of the
            %   derivatives D(l,:) of order at most M <= 4 of a product
            %   f*g of two functions whose derivatives are stored as
            %   columns in the order of the rows of
            %   multiIndices(size(D,2),M): the rows TERMS.l of D, the
            %   columns TERMS.f and TERMS.g of the derivatives of f and
            %   g, the multinomial coefficients TERMS.coeff, and the
            %   number TERMS.n of rows of D.
            dim = size(D,2);
            Dall = BarycentricJacobi.multiIndices(dim,m);
            key = (m+1).^(0:dim-1)';
            row4D = zeros((m+1)^dim,1);
            row4D(Dall*key+1) = 1:size(Dall,1);
            fact = factorial(0:m);
            terms = struct('l',[],'f',[],'g',[],'coeff',[],'n',size(D,1));
            for l = 1 : size(D,1)
                for b = find(all(Dall <= D(l,:),2))'
                    gamma = D(l,:) - Dall(b,:);
                    terms.l(end+1) = l;
                    terms.f(end+1) = b;
                    terms.g(end+1) = row4D(gamma*key+1);
                    terms.coeff(end+1) = prod(fact(D(l,:)+1)./(fact(Dall(b,:)+1).*fact(gamma+1)));
                end
            end
        end

        function val = product(f,G,terms)
            %PRODUCT - Derivatives of products of two factors
            %   VAL = PRODUCT(F,G,TERMS) evaluates the derivatives of the
            %   products f*g_k by the Leibniz rule TERMS of leibnizTerms,
            %   where the vector F and the rows of G hold the derivatives
            %   of f and g_k in the order of multiIndices. VAL has one
            %   row per g_k and one column per row of the multi-indices
            %   D of leibnizTerms.
            val = (terms.coeff.*reshape(f(terms.f),1,[]).*G(:,terms.g))*double(terms.l' == 1:terms.n);
        end

        function [w,curlw] = whitney(lambda,gradLambda,a,b)
            %WHITNEY - Whitney edge function and its curl
            %   W = WHITNEY(LAMBDA,GRADLAMBDA,A,B) returns the Whitney
            %   function W = lambda_A grad(lambda_B) - lambda_B
            %   grad(lambda_A) (1-by-d) of the edge [A B] at the
            %   barycentric coordinates LAMBDA of one point with the
            %   gradients GRADLAMBDA from coordinates.
            %
            %   [W,CURLW] = WHITNEY(LAMBDA,GRADLAMBDA,A,B) also returns
            %   the constant curl CURLW = 2 grad(lambda_A) x
            %   grad(lambda_B) (1-by-3) in 3d or the rot
            %   CURLW = 2 det([grad(lambda_A); grad(lambda_B)]) in 2d.
            w = lambda(a)*gradLambda(b,:) - lambda(b)*gradLambda(a,:);
            if size(gradLambda,2) == 3
                curlw = 2*BarycentricJacobi.cross(gradLambda(a,:),gradLambda(b,:));
            else
                curlw = 2*det(gradLambda([a b],:));
            end
        end

        function [phi,divphi] = whitneyFace(lambda,gradLambda,a,b,c)
            %whitneyFace - Whitney face function and its divergence
            %   PHI = whitneyFace(LAMBDA,GRADLAMBDA,A,B,C) returns the
            %   Whitney function PHI = lambda_A grad(lambda_B) x
            %   grad(lambda_C) + lambda_B grad(lambda_C) x grad(lambda_A)
            %   + lambda_C grad(lambda_A) x grad(lambda_B) (1-by-3) of
            %   the face [A B C] of the reference tetrahedron at the
            %   barycentric coordinates LAMBDA of one point with the
            %   gradients GRADLAMBDA from coordinates.
            %
            %   [PHI,DIVPHI] = whitneyFace(LAMBDA,GRADLAMBDA,A,B,C) also
            %   returns the constant divergence DIVPHI =
            %   3 det([grad(lambda_A); grad(lambda_B); grad(lambda_C)]).
            g = gradLambda([a b c],:);
            phi = lambda(a)*BarycentricJacobi.cross(g(2,:),g(3,:)) ...
                + lambda(b)*BarycentricJacobi.cross(g(3,:),g(1,:)) ...
                + lambda(c)*BarycentricJacobi.cross(g(1,:),g(2,:));
            divphi = 3*det(g);
        end

        function val = timesField(F,w,dw,op)
            %timesField - Products of functions with a vector field
            %   VAL = timesField(F,W,DW,OP) returns the products f_k W of
            %   the functions f_k with the vector field W for OP = '',
            %   their curls grad(f_k) x W + f_k DW in 3d or their rots in
            %   2d for OP = 'curl', or their divergences grad(f_k).W +
            %   f_k DW for OP = 'div'. The rows of F hold the values and
            %   gradients [f_k, grad(f_k)] (d+1 columns), W is 1-by-d,
            %   and DW is the curl, rot, or divergence of W at the same
            %   point. VAL has one row per function f_k: N-by-d for
            %   OP = '' and for the curls in 3d, and N-by-1 for the rots
            %   and divergences, where N = size(F,1).
            if isempty(op)
                val = F(:,1).*w;
            elseif strcmp(op,'div')
                val = F(:,2:end)*w(:) + F(:,1)*dw;
            elseif numel(w) == 3
                val = BarycentricJacobi.cross(F(:,2:4),w) + F(:,1).*dw;
            else
                val = F(:,2)*w(2) - F(:,3)*w(1) + F(:,1)*dw;
            end
        end

        function c = cross(a,b)
            %CROSS - Cross products of rows
            %   C = CROSS(A,B) returns the cross products of the rows of
            %   A and B with three columns, where a single row is used
            %   for all rows of the other argument.
            c = [a(:,2).*b(:,3)-a(:,3).*b(:,2), a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
                a(:,1).*b(:,2)-a(:,2).*b(:,1)];
        end

        function D = multiIndices(dim,m)
            %multiIndices - Multi-indices of all derivatives up to an order
            %   D = multiIndices(DIM,M) returns the multi-indices of the
            %   partial derivatives of the orders 0,...,M <= 4 in DIM = 2
            %   or 3 variables, one row per derivative, ordered by their
            %   order and within an order as in Monomials.multiIndices:
            %   the value, the gradient, the second derivatives, and so
            %   on.
            D = zeros(1,dim);
            for k = 1 : m
                D = [D; Monomials.multiIndices(dim,k)]; %#ok<AGROW>
            end
        end

        function cols = columns(dim,order)
            %COLUMNS - Columns of the derivatives of one order
            %   COLS = COLUMNS(DIM,ORDER) returns the columns of the
            %   derivatives of order ORDER, i.e. the indices of their
            %   rows in multiIndices(DIM,M) for any M >= ORDER.
            first = 0;
            if order > 0
                first = nchoosek(order-1+dim,dim);
            end
            cols = first+1 : nchoosek(order+dim,dim);
        end

    end

end
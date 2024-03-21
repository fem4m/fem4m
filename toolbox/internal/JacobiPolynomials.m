classdef JacobiPolynomials
    %JacobiPolynomials - Jacobi polynomials and their integrals
    %   JacobiPolynomials provides static methods that evaluate the
    %   Jacobi polynomials P_j^(ALPHA,BETA) on [-1,1] of degree
    %   j = 0,...,N, their integrated and twice integrated versions
    %   (from -1, with BETA = 0), the scaled versions S^j P_j(2X/S-1)
    %   of the Jacobi and integrated Jacobi polynomials, or their
    %   derivatives, returning one row per degree and one column per
    %   point.
    %
    %   JacobiPolynomials functions:
    %       twiceIntegratedPoly  - Twice integrated Jacobi polynomials
    %       integratedPoly       - Integrated Jacobi polynomials
    %       poly                 - Jacobi polynomials
    %       scaledPoly           - Scaled Jacobi polynomials
    %       scaledIntegratedPoly - Scaled integrated Jacobi polynomials
    %
    %   See also BarycentricJacobi, Monomials, HPk, DGDk, AWk

    % Copyright (c) 2026 by the fem4m authors

    methods (Static)

        function P = twiceIntegratedPoly(x,n,alpha,k)
            %twiceIntegratedPoly - Twice integrated Jacobi polynomials
            %   P = twiceIntegratedPoly(X,N,ALPHA,K) returns the K-th
            %   derivatives of the twice integrated Jacobi polynomials
            %   M_j(X) = int_{-1}^X L_{j-1}(t) dt, j = 1,...,N, and
            %   M_0 = 1, with the integrated Jacobi polynomials L_j of
            %   integratedPoly, at the points of the column vector X,
            %   one row per degree j and one column per point.
            P = zeros(n+1,numel(x));
            J = JacobiPolynomials.integratedPoly(x,n,alpha,k);
            if k == 0;P(1,:) = 1;end
            if n >= 1
                if k == 0;P(2,:) = x+1;elseif k == 1; P(2,:) = 1;end
            end
            if n >= 2
                if k == 0;P(3,:) = .5*x.^2+x-.5;elseif k == 1; P(3,:) = x+1;elseif k == 2; P(3,:) = 1;end
            end
            for j = 2:n-1
                a1 = (2*j+2*alpha)/(2*j+alpha-1)/(2*j+alpha);
                a2 = 2*alpha/(2*j+alpha-2)/(2*j+alpha);
                a3 = (2*j-2)/(2*j+alpha-1)/(2*j+alpha-2);
                P(j+2,:) = a1*J(j+2,:) + a2*J(j+1,:) - a3*J(j,:);
            end
        end

        function P = integratedPoly(x,n,alpha,k)
            %integratedPoly - Integrated Jacobi polynomials
            %   P = integratedPoly(X,N,ALPHA,K) returns the K-th
            %   derivatives of the integrated Jacobi polynomials
            %   L_j(X) = int_{-1}^X P_{j-1}^(ALPHA,0)(t) dt, j = 1,...,N,
            %   and L_0 = 1 at the points of the column vector X, one
            %   row per degree j and one column per point.
            P = zeros(n+1,numel(x));
            J = JacobiPolynomials.poly(x,n,alpha,0,k);
            if k == 0;P(1,:) = 1;end
            if n >= 1
                if k == 0;P(2,:) = x+1;elseif k == 1; P(2,:) = 1;end
            end
            for j = 2:n
                a1 = (2*j+2*alpha)/(2*j+alpha-1)/(2*j+alpha);
                a2 = 2*alpha/(2*j+alpha-2)/(2*j+alpha);
                a3 = (2*j-2)/(2*j+alpha-1)/(2*j+alpha-2);
                P(j+1,:) = a1*J(j+1,:) + a2*J(j,:) - a3*J(j-1,:);
            end
        end

        function P = poly(x,n,alpha,beta,k)
            %POLY - Jacobi polynomials
            %   P = POLY(X,N,ALPHA,BETA,K) returns the K-th derivatives
            %   of the Jacobi polynomials P_j^(ALPHA,BETA), j = 0,...,N,
            %   normalized by P_j(1) = nchoosek(j+ALPHA,j), at the points
            %   of the column vector X, one row per degree j and one
            %   column per point; the rows of degree j < K are zero.
            P = zeros(n+1,numel(x));
            % formula with k-th derivative shift in n, alpha and beta
            if k+1 <= n+1
                P(k+1,:) = 1;
            end
            if n-k > 0
                P(k+2,:) = 1/2*(alpha-beta+(alpha+beta+2*(1+k))*x');
                for j = 1 : n-k-1
                    a1 = 2*(j+1)*(j+alpha+beta+1+2*k)*(2*j+alpha+beta+2*k);
                    a2 = (2*j+alpha+beta+1+2*k)*((alpha+k)^2-(beta+k)^2);
                    a3 = (2*j+alpha+beta+2*k)*(2*j+alpha+beta+1+2*k)*(2*j+alpha+beta+2+2*k);
                    a4 = 2*(j+alpha+k)*(j+beta+k)*(2*j+alpha+beta+2+2*k);
                    P(k+j+2,:) = ( (a2+a3*x').*P(k+j+1,:) - a4*P(k+j,:) )/a1;
                end
            end
            % derivative scaling
            for j = k : n
                P(j+1,:) = P(j+1,:)*gamma(alpha+beta+j+1+k)*2^(-k)/gamma(alpha+beta+j+1);
            end
        end

        function Q = scaledPoly(x,s,n,alpha,beta,kx,ks)
            %scaledPoly - Scaled Jacobi polynomials
            %   Q = scaledPoly(X,S,N,ALPHA,BETA,KX,KS) returns the
            %   derivatives, KX times in X and KS times in S, of the
            %   scaled Jacobi polynomials Q_j(X,S) = S^j P_j(2X/S-1),
            %   j = 0,...,N, with the Jacobi polynomials P_j of poly with
            %   parameters ALPHA and BETA, at the points (X,S), one row
            %   per degree j and one column per point; the rows of degree
            %   j < KX+KS are zero. X and S are vectors of equal length,
            %   or S is a scalar. The Q_j are polynomials in X and S, also
            %   for S = 0.
            x = x(:)';
            s = s(:)';
            k = kx+ks;
            Q = zeros(n+1,numel(x));
            % d/dx Q_j^(a,b) = (j+a+b+1) Q_{j-1}^(a+1,b+1) and
            % d/ds Q_j^(a,b) = -(j+b) Q_{j-1}^(a+1,b) shift a and b
            a = alpha+k;
            b = beta+kx;
            if k+1 <= n+1
                Q(k+1,:) = 1;
            end
            if n-k > 0
                Q(k+2,:) = (a+b+2)*x-(b+1)*s;
                for j = 1 : n-k-1
                    a1 = 2*(j+1)*(j+a+b+1)*(2*j+a+b);
                    a2 = (2*j+a+b+1)*(a^2-b^2);
                    a3 = (2*j+a+b)*(2*j+a+b+1)*(2*j+a+b+2);
                    a4 = 2*(j+a)*(j+b)*(2*j+a+b+2);
                    Q(k+j+2,:) = ( ((a2-a3)*s+2*a3*x).*Q(k+j+1,:) - a4*s.^2.*Q(k+j,:) )/a1;
                end
            end
            % derivative scaling
            for j = k : n
                Q(j+1,:) = Q(j+1,:)*(-1)^ks*prod(j+beta-(0:ks-1))*prod(j+alpha+beta+1+(0:kx-1));
            end
        end

        function L = scaledIntegratedPoly(x,s,n,alpha,kx,ks)
            %scaledIntegratedPoly - Scaled integrated Jacobi polynomials
            %   L = scaledIntegratedPoly(X,S,N,ALPHA,KX,KS) returns the
            %   derivatives, KX times in X and KS times in S, of the
            %   scaled integrated Jacobi polynomials
            %   L_j(X,S) = S^j L_j(2X/S-1), j = 0,...,N, with the
            %   integrated Jacobi polynomials L_j of integratedPoly
            %   (from -1, with BETA = 0), at the points (X,S), one row
            %   per degree j and one column per point. X and S are
            %   vectors of equal length, or S is a scalar. The L_j are
            %   polynomials in X and S, also for S = 0.
            x = x(:)';
            s = s(:)';
            L = zeros(n+1,numel(x));
            if kx > 0
                % d/dx L_j = 2 Q_{j-1}^(alpha,0)
                L(2:end,:) = 2*JacobiPolynomials.scaledPoly(x,s,n-1,alpha,0,kx-1,ks);
            elseif ks > 0
                % d/ds L_j = -2(j-1)/(j+alpha-1) Q_{j-1}^(alpha,-1)
                Q = JacobiPolynomials.scaledPoly(x,s,n-1,alpha,-1,0,ks-1);
                j = (2:n)';
                L(3:end,:) = -2*(j-1)./(j+alpha-1).*Q(2:end,:);
            else
                % L_j = a1*Q_j + a2*s*Q_{j-1} - a3*s^2*Q_{j-2}
                Q = JacobiPolynomials.scaledPoly(x,s,n,alpha,0,0,0);
                L(1,:) = 1;
                if n >= 1
                    L(2,:) = 2*x;
                end
                for j = 2:n
                    a1 = (2*j+2*alpha)/(2*j+alpha-1)/(2*j+alpha);
                    a2 = 2*alpha/(2*j+alpha-2)/(2*j+alpha);
                    a3 = (2*j-2)/(2*j+alpha-1)/(2*j+alpha-2);
                    L(j+1,:) = a1*Q(j+1,:) + a2*s.*Q(j,:) - a3*s.^2.*Q(j-1,:);
                end
            end
        end
    end
end
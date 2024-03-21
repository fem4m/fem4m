classdef problemScalarEllipticEigenvalueSquareTwist < ProblemEigenvalue
    %problemScalarEllipticEigenvalueSquareTwist - Twisted diffusion on the unit square
    %   OBJ = problemScalarEllipticEigenvalueSquareTwist() supplies the
    %   Dirichlet eigenvalue problem -div(alpha*grad(u)) = lambda*u on the
    %   unit square with a smooth, full and anisotropic diffusion matrix
    %   alpha and the analytic eigenvalues pi^2*(j^2+k^2), j,k >= 1. The
    %   modes M = 1,2,... number the eigenvalues in increasing order,
    %   repeated according to their multiplicity.
    %
    %   OBJ = problemScalarEllipticEigenvalueSquareTwist(THETA0) sets the
    %   maximal twist angle THETA0 (default pi/2).
    %
    %   The area preserving twist Phi(x) = c + Rot(g(r))*(x-c), with the
    %   rotation Rot(g) by the angle g, r = |x-c|, c = (1/2,1/2) and
    %   g(r) = THETA0*(1-r^2/R^2)^4 for r < R = 1/2 (g = 0 otherwise), maps
    %   the unit square onto itself with det(DPhi) = 1. With
    %   alpha = inv(DPhi)*inv(DPhi)' the functions u = v(Phi(x)) are the
    %   eigenfunctions, where v(y) = 2*sin(j*pi*y_1)*sin(k*pi*y_2) are the
    %   Dirichlet Laplace eigenfunctions of the unit square; alpha is C^2,
    %   symmetric positive definite with det(alpha) = 1 and equals the
    %   identity outside the disk r < R.
    %
    %   problemScalarEllipticEigenvalueSquareTwist functions:
    %       eigenvalue  - Analytic eigenvalue(s) of mode(s) M
    %       eigenspace  - Index pairs [j k] of the eigenspace of mode M
    %       cluster     - Modes with the same eigenvalue as mode M
    %       alpha       - Diffusion matrix
    %       dalpha      - Derivative of the diffusion matrix
    %       Phi         - Twist map
    %       DPhi        - Jacobian of the twist map
    %       exactU      - L^2-normalized eigenfunction of pair (J,K)
    %       exactGradU  - Gradient of the eigenfunction of pair (J,K)
    %
    %   problemScalarEllipticEigenvalueSquareTwist properties:
    %       geom        - Unit square geometry
    %       theta0      - Maximal twist angle at the center
    %       radius      - Radius of the twisted disk
    %       center      - Center of the twisted disk
    %
    %   See also ProblemEigenvalue, problemScalarEllipticEigenvalueSquare

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        theta0 = pi/2
        radius = 1/2
        center = [1/2 1/2]
    end

    methods

        function obj = problemScalarEllipticEigenvalueSquareTwist(theta0)
            if nargin > 0
                obj.theta0 = theta0;
            end
        end

        function val = eigenvalue(obj,m)
            val = pi^2*obj.squareSums(m);
        end

        function val = eigenspace(obj,m)
            %EIGENSPACE - Index pairs [j k] of the eigenspace of mode M
            %   VAL = EIGENSPACE(OBJ,M) returns the pairs [j k] with
            %   j^2+k^2 = eigenvalue(M)/pi^2, one row for each
            %   L^2-orthonormal eigenfunction exactU(x,j,k). M must be a
            %   scalar.
            n = obj.squareSums(m);
            [j,k] = ndgrid(1:floor(sqrt(n)));
            ind = j(:).^2 + k(:).^2 == n;
            val = [j(ind) k(ind)];
        end

        function val = cluster(obj,m)
            %CLUSTER - Modes with the same eigenvalue as mode M
            %   VAL = CLUSTER(OBJ,M) returns the modes (a row vector) whose
            %   eigenvalue equals eigenvalue(M). M must be a scalar.
            n = obj.squareSums(m);
            val = find(obj.squareSums(1:n) == n)';
        end

        function val = alpha(obj,x)
            [xi,eta,kappa] = obj.twist(x);
            r2 = xi.^2 + eta.^2;
            val = zeros(size(x,1),2,2);
            val(:,1,1) = 1 + 2*kappa.*xi.*eta + kappa.^2.*r2.*eta.^2;
            val(:,1,2) = -kappa.*(xi.^2-eta.^2) - kappa.^2.*r2.*xi.*eta;
            val(:,2,1) = val(:,1,2);
            val(:,2,2) = 1 - 2*kappa.*xi.*eta + kappa.^2.*r2.*xi.^2;
        end

        function val = dalpha(obj,x)
            %DALPHA - Derivative of the diffusion matrix
            %   VAL = DALPHA(OBJ,X) returns the derivatives at the points X
            %   (N-by-2, one row per point) as an N-by-2-by-2-by-2 array,
            %   VAL(:,i,j,k) = d_k alpha_ij.
            [xi,eta,kappa,dkappa] = obj.twist(x);
            r2 = xi.^2 + eta.^2;
            val = zeros(size(x,1),2,2,2);
            for k = 1:2
                dk = dkappa(:,k);
                if k == 1
                    dxieta = eta; dxi2eta2 = 2*xi;
                    dr2eta2 = 2*xi.*eta.^2;
                    dr2xieta = 2*xi.^2.*eta + r2.*eta;
                    dr2xi2 = 2*xi.^3 + 2*xi.*r2;
                else
                    dxieta = xi; dxi2eta2 = -2*eta;
                    dr2eta2 = 2*eta.^3 + 2*eta.*r2;
                    dr2xieta = 2*xi.*eta.^2 + r2.*xi;
                    dr2xi2 = 2*eta.*xi.^2;
                end
                val(:,1,1,k) = 2*dk.*xi.*eta + 2*kappa.*dxieta ...
                    + 2*kappa.*dk.*r2.*eta.^2 + kappa.^2.*dr2eta2;
                val(:,1,2,k) = -dk.*(xi.^2-eta.^2) - kappa.*dxi2eta2 ...
                    - 2*kappa.*dk.*r2.*xi.*eta - kappa.^2.*dr2xieta;
                val(:,2,1,k) = val(:,1,2,k);
                val(:,2,2,k) = -2*dk.*xi.*eta - 2*kappa.*dxieta ...
                    + 2*kappa.*dk.*r2.*xi.^2 + kappa.^2.*dr2xi2;
            end
        end

        function val = Phi(obj,x)
            %PHI - Twist map
            %   VAL = PHI(OBJ,X) returns the images Phi(X) (N-by-2) of the
            %   points X (N-by-2, one row per point).
            [xi,eta,~,~,g] = obj.twist(x);
            val = obj.center + [cos(g).*xi - sin(g).*eta, sin(g).*xi + cos(g).*eta];
        end

        function val = DPhi(obj,x)
            %DPhi - Jacobian of the twist map
            %   VAL = DPhi(OBJ,X) returns the Jacobian at the points X
            %   (N-by-2, one row per point) as an N-by-2-by-2 array,
            %   VAL(:,i,k) = d_k Phi_i. It is
            %   DPhi = Rot(g)*(I + kappa*J*(x-c)*(x-c)^T) with kappa = g'(r)/r,
            %   the identity I and J = [0 -1; 1 0].
            [xi,eta,kappa,~,g] = obj.twist(x);
            M11 = 1 - kappa.*xi.*eta;
            M12 = -kappa.*eta.^2;
            M21 = kappa.*xi.^2;
            M22 = 1 + kappa.*xi.*eta;
            c = cos(g);
            s = sin(g);
            val = zeros(size(x,1),2,2);
            val(:,1,1) = c.*M11 - s.*M21;
            val(:,1,2) = c.*M12 - s.*M22;
            val(:,2,1) = s.*M11 + c.*M21;
            val(:,2,2) = s.*M12 + c.*M22;
        end

        function val = exactU(obj,x,j,k)
            %exactU - L^2-normalized eigenfunction of pair (J,K)
            %   VAL = exactU(OBJ,X) returns the N-by-1 values of
            %   u = 2*sin(pi*Phi_1)*sin(pi*Phi_2), the pair (J,K) = (1,1), at the
            %   points X (N-by-2, one row per point).
            %
            %   VAL = exactU(OBJ,X,J,K) uses the pair (J,K),
            %   u = 2*sin(J*pi*Phi_1)*sin(K*pi*Phi_2).
            if nargin < 3, j = 1; k = 1; end
            y = obj.Phi(x);
            val = 2*sin(j*pi*y(:,1)).*sin(k*pi*y(:,2));
        end

        function val = exactGradU(obj,x,j,k)
            %exactGradU - Gradient of the eigenfunction of pair (J,K)
            %   VAL = exactGradU(OBJ,X) returns the N-by-2 gradients
            %   DPhi'*grad(v)(Phi) of u = v(Phi) at the points X (N-by-2, one row
            %   per point) for v = 2*sin(pi*y_1)*sin(pi*y_2), the pair (1,1).
            %
            %   VAL = exactGradU(OBJ,X,J,K) uses the pair (J,K),
            %   v = 2*sin(J*pi*y_1)*sin(K*pi*y_2).
            if nargin < 3, j = 1; k = 1; end
            y = obj.Phi(x);
            D = obj.DPhi(x);
            gradV = 2*pi*[j*cos(j*pi*y(:,1)).*sin(k*pi*y(:,2)), ...
                k*sin(j*pi*y(:,1)).*cos(k*pi*y(:,2))];
            val = [gradV(:,1).*D(:,1,1) + gradV(:,2).*D(:,2,1), ...
                gradV(:,1).*D(:,1,2) + gradV(:,2).*D(:,2,2)];
        end

    end

    methods (Access = private)

        function [xi,eta,kappa,dkappa,g] = twist(obj,x)
            % local coordinates x-c, kappa = g'(r)/r, its gradient and the angle g(r)
            xi = x(:,1) - obj.center(1);
            eta = x(:,2) - obj.center(2);
            R2 = obj.radius^2;
            w = max(1 - (xi.^2+eta.^2)/R2,0);
            kappa = -8*obj.theta0/R2*w.^3;
            dkappa = 48*obj.theta0/R2^2*w.^2.*[xi eta];
            g = obj.theta0*w.^4;
        end

    end

    methods (Static, Access = private)

        function val = squareSums(m)
            % m-th smallest value of j^2+k^2 with j,k >= 1 (with multiplicity)
            j = 1:max(m);
            k = (1:max(m))';
            val = sort(reshape(j.^2+k.^2,[],1));
            val = val(m);
        end

    end

end
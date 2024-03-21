classdef B2k < DGDk
    %B2k - Squared cubic bubble times Dubiner polynomials of degree k-6
    %   OBJ = B2k(DEGREE,MESH) creates the element bubble space
    %   of polynomial degree DEGREE >= 6 (a scalar or one value per
    %   element) on the triangular mesh MESH whose element basis
    %   functions are the Dubiner polynomials of DGDk of total degree
    %   DEGREE-6 multiplied by the squared cubic bubble
    %   b_T^2 = (x*y*(1-x-y))^2 of the reference triangle. The basis
    %   functions and their gradients vanish on the boundary of their
    %   element.
    %
    %   B2k functions:
    %       getElementBasis     - Evaluate the bubble element basis
    %       getGradElementBasis - Evaluate the gradient of the bubbles
    %       getD2ElementBasis   - Evaluate the Hessian of the bubbles
    %       getD3ElementBasis   - Evaluate third derivatives of bubbles
    %       getD4ElementBasis   - Evaluate fourth derivatives of bubbles
    %
    %   See also Bk, CurlBk, DGDk, Pk

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = B2k(degree,mesh)
            %B2k - Create a squared cubic bubble element space
            if any(degree(:) < 6)
                error('fem4m:B2k:invalidDegree','B2k degree must be at least 6');
            end
            obj = obj@DGDk(degree-6,mesh);
            obj.degree4e = obj.degree4e+6;
            obj.degree = obj.degree+6;
            obj.degree4DoF = obj.degree4DoF+6;
        end
    end

    methods (Static)
        %% Elem basis helper functions
        function val = getElementBasis(p,x)
            %getElementBasis - Evaluate the bubble element basis
            %   VAL = getElementBasis(P,X) multiplies the Dubiner
            %   basis of total degree P-6 with the squared cubic
            %   bubble at the reference point X (1-by-2). VAL is
            %   1-by-(P-5)(P-4)/2.
            basis = DGDk.getElementBasis(p-6,x);
            y = x(:,2)';
            x = x(:,1)';
            val = ((1-x-y).*x.*y).^2.*basis;
        end

        function val = getGradElementBasis(p,x)
            %getGradElementBasis - Evaluate the gradient of the bubbles
            %   VAL = getGradElementBasis(P,X) evaluates the
            %   gradients of the squared cubic bubble times the
            %   Dubiner basis of total degree P-6 at the reference
            %   point X (1-by-2). VAL is (P-5)(P-4)/2-by-2.
            basis = DGDk.getElementBasis(p-6,x)';
            grad = DGDk.getGradElementBasis(p-6,x);
            y = x(:,2)';
            x = x(:,1)';
            val(:,1,:) = x.^2.*y.^2.*basis.*(x.*2.0+y.*2.0-2.0)+x.*y.^2.*basis.*(x+y-1.0).^2.*2.0+x.^2.*y.^2.*grad(:,1,:).*(x+y-1.0).^2;
            val(:,2,:) = x.^2.*y.^2.*basis.*(x.*2.0+y.*2.0-2.0)+x.^2.*y.*basis.*(x+y-1.0).^2.*2.0+x.^2.*y.^2.*grad(:,2,:).*(x+y-1.0).^2;
        end

        function val = getD2ElementBasis(p,x)
            %getD2ElementBasis - Evaluate the Hessian of the bubbles
            %   VAL = getD2ElementBasis(P,X) evaluates the Hessians
            %   [dxx,dyy,dxy] of the squared cubic bubble times the
            %   Dubiner basis of total degree P-6 at the reference
            %   point X (1-by-2). VAL is (P-5)(P-4)/2-by-3.
            basis = DGDk.getElementBasis(p-6,x)';
            grad = DGDk.getGradElementBasis(p-6,x);
            d2 = DGDk.getD2ElementBasis(p-6,x);
            y = x(:,2)';
            x = x(:,1)';
            val(:,1,:) = y.^2.*basis.*(x+y-1.0).^2.*2.0+x.^2.*y.^2.*basis.*2.0+x.*y.^2.*basis.*(x.*2.0+y.*2.0-2.0).*4.0+x.^2.*y.^2.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*2.0+x.^2.*y.^2.*(x+y-1.0).^2.*d2(:,1,:)+x.*y.^2.*grad(:,1,:).*(x+y-1.0).^2.*4.0;
            val(:,2,:) = x.^2.*basis.*(x+y-1.0).^2.*2.0+x.^2.*y.^2.*basis.*2.0+x.^2.*y.*basis.*(x.*2.0+y.*2.0-2.0).*4.0+x.^2.*y.^2.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*2.0+x.^2.*y.^2.*(x+y-1.0).^2.*d2(:,2,:)+x.^2.*y.*grad(:,2,:).*(x+y-1.0).^2.*4.0;
            val(:,3,:) = x.^2.*y.^2.*basis.*2.0+x.*y.^2.*basis.*(x.*2.0+y.*2.0-2.0).*2.0+x.^2.*y.*basis.*(x.*2.0+y.*2.0-2.0).*2.0+x.^2.*y.^2.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0)+x.^2.*y.^2.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0)+x.*y.*basis.*(x+y-1.0).^2.*4.0+x.^2.*y.^2.*(x+y-1.0).^2.*d2(:,3,:)+x.^2.*y.*grad(:,1,:).*(x+y-1.0).^2.*2.0+x.*y.^2.*grad(:,2,:).*(x+y-1.0).^2.*2.0;
        end

        function val = getD3ElementBasis(p,x)
            %getD3ElementBasis - Evaluate third derivatives of bubbles
            %   VAL = getD3ElementBasis(P,X) evaluates the third
            %   derivatives [u_xxx,u_yyy,u_xxy,u_xyy] of the squared
            %   cubic bubble times the Dubiner basis of total degree
            %   P-6 at the reference point X (1-by-2). VAL is
            %   (P-5)(P-4)/2-by-4.
            basis = DGDk.getElementBasis(p-6,x)';
            grad = DGDk.getGradElementBasis(p-6,x);
            d2 = DGDk.getD2ElementBasis(p-6,x);
            d3 = DGDk.getD3ElementBasis(p-6,x);
            y = x(:,2)';
            x = x(:,1)';
            val(:,1,:) = x.*y.^2.*basis.*1.2e+1+y.^2.*grad(:,1,:).*(x+y-1.0).^2.*6.0+x.^2.*y.^2.*grad(:,1,:).*6.0+y.^2.*basis.*(x.*2.0+y.*2.0-2.0).*6.0+x.*y.^2.*(x+y-1.0).^2.*d2(:,1,:).*6.0+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,1,:).*3.0+x.^2.*y.^2.*(x+y-1.0).^2.*d3(:,1,:)+x.*y.^2.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*1.2e+1;
            val(:,2,:) = x.^2.*y.*basis.*1.2e+1+x.^2.*grad(:,2,:).*(x+y-1.0).^2.*6.0+x.^2.*y.^2.*grad(:,2,:).*6.0+x.^2.*basis.*(x.*2.0+y.*2.0-2.0).*6.0+x.^2.*y.*(x+y-1.0).^2.*d2(:,2,:).*6.0+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,2,:).*3.0+x.^2.*y.^2.*(x+y-1.0).^2.*d3(:,2,:)+x.^2.*y.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*1.2e+1;
            val(:,3,:) = x.*y.^2.*basis.*8.0+x.^2.*y.*basis.*4.0+y.^2.*grad(:,2,:).*(x+y-1.0).^2.*2.0+x.^2.*y.^2.*grad(:,1,:).*4.0+x.^2.*y.^2.*grad(:,2,:).*2.0+y.^2.*basis.*(x.*2.0+y.*2.0-2.0).*2.0+y.*basis.*(x+y-1.0).^2.*4.0+x.^2.*y.*(x+y-1.0).^2.*d2(:,1,:).*2.0+x.*y.^2.*(x+y-1.0).^2.*d2(:,3,:).*4.0+x.*y.*grad(:,1,:).*(x+y-1.0).^2.*8.0+x.*y.*basis.*(x.*2.0+y.*2.0-2.0).*8.0+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,1,:)+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,3,:).*2.0+x.^2.*y.^2.*(x+y-1.0).^2.*d3(:,3,:)+x.*y.^2.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*4.0+x.^2.*y.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*4.0+x.*y.^2.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*4.0;
            val(:,4,:) = x.*y.^2.*basis.*4.0+x.^2.*y.*basis.*8.0+x.^2.*grad(:,1,:).*(x+y-1.0).^2.*2.0+x.^2.*y.^2.*grad(:,1,:).*2.0+x.^2.*y.^2.*grad(:,2,:).*4.0+x.^2.*basis.*(x.*2.0+y.*2.0-2.0).*2.0+x.*basis.*(x+y-1.0).^2.*4.0+x.^2.*y.*(x+y-1.0).^2.*d2(:,3,:).*4.0+x.*y.^2.*(x+y-1.0).^2.*d2(:,2,:).*2.0+x.*y.*grad(:,2,:).*(x+y-1.0).^2.*8.0+x.*y.*basis.*(x.*2.0+y.*2.0-2.0).*8.0+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,3,:).*2.0+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,2,:)+x.^2.*y.^2.*(x+y-1.0).^2.*d3(:,4,:)+x.^2.*y.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*4.0+x.*y.^2.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*4.0+x.^2.*y.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*4.0;
        end

        function val = getD4ElementBasis(p,x)
            %getD4ElementBasis - Evaluate fourth derivatives of bubbles
            %   VAL = getD4ElementBasis(P,X) evaluates the fourth
            %   derivatives [u_xxxx,u_yyyy,u_xxyy,u_xxxy,u_xyyy] of
            %   the squared cubic bubble times the Dubiner basis of
            %   total degree P-6 at the reference point X (1-by-2).
            %   VAL is (P-5)(P-4)/2-by-5.
            basis = DGDk.getElementBasis(p-6,x)';
            grad = DGDk.getGradElementBasis(p-6,x);
            d2 = DGDk.getD2ElementBasis(p-6,x);
            d3 = DGDk.getD3ElementBasis(p-6,x);
            d4 = DGDk.getD4ElementBasis(p-6,x);
            y = x(:,2)';
            x = x(:,1)';
            val(:,1,:) = y.^2.*basis.*2.4e+1+y.^2.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*2.4e+1+y.^2.*(x+y-1.0).^2.*d2(:,1,:).*1.2e+1+x.^2.*y.^2.*d2(:,1,:).*1.2e+1+x.*y.^2.*grad(:,1,:).*4.8e+1+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d3(:,1,:).*4.0+x.^2.*y.^2.*(x+y-1.0).^2.*d4(:,1,:)+x.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,1,:).*2.4e+1+x.*y.^2.*(x+y-1.0).^2.*d3(:,1,:).*8.0;
            val(:,2,:) = x.^2.*basis.*2.4e+1+x.^2.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*2.4e+1+x.^2.*(x+y-1.0).^2.*d2(:,2,:).*1.2e+1+x.^2.*y.^2.*d2(:,2,:).*1.2e+1+x.^2.*y.*grad(:,2,:).*4.8e+1+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d3(:,2,:).*4.0+x.^2.*y.^2.*(x+y-1.0).^2.*d4(:,2,:)+x.^2.*y.*(x.*2.0+y.*2.0-2.0).*d2(:,2,:).*2.4e+1+x.^2.*y.*(x+y-1.0).^2.*d3(:,2,:).*8.0;
            val(:,3,:) = basis.*(x+y-1.0).^2.*4.0+x.^2.*basis.*4.0+y.^2.*basis.*4.0+x.*basis.*(x.*2.0+y.*2.0-2.0).*8.0+y.*basis.*(x.*2.0+y.*2.0-2.0).*8.0+x.^2.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*4.0+y.^2.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*4.0+x.^2.*(x+y-1.0).^2.*d2(:,1,:).*2.0+y.^2.*(x+y-1.0).^2.*d2(:,2,:).*2.0+x.*y.*basis.*3.2e+1+x.*grad(:,1,:).*(x+y-1.0).^2.*8.0+x.^2.*y.^2.*d2(:,1,:).*2.0+x.^2.*y.^2.*d2(:,3,:).*8.0+y.*grad(:,2,:).*(x+y-1.0).^2.*8.0+x.^2.*y.^2.*d2(:,2,:).*2.0+x.*y.^2.*grad(:,1,:).*8.0+x.^2.*y.*grad(:,1,:).*1.6e+1+x.*y.^2.*grad(:,2,:).*1.6e+1+x.^2.*y.*grad(:,2,:).*8.0+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d3(:,3,:).*2.0+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d3(:,4,:).*2.0+x.^2.*y.^2.*(x+y-1.0).^2.*d4(:,3,:)+x.^2.*y.*(x.*2.0+y.*2.0-2.0).*d2(:,1,:).*4.0+x.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,3,:).*8.0+x.^2.*y.*(x.*2.0+y.*2.0-2.0).*d2(:,3,:).*8.0+x.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,2,:).*4.0+x.^2.*y.*(x+y-1.0).^2.*d3(:,3,:).*4.0+x.*y.^2.*(x+y-1.0).^2.*d3(:,4,:).*4.0+x.*y.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*1.6e+1+x.*y.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*1.6e+1+x.*y.*(x+y-1.0).^2.*d2(:,3,:).*1.6e+1;
            val(:,4,:) = y.^2.*basis.*1.2e+1+y.*basis.*(x.*2.0+y.*2.0-2.0).*1.2e+1+y.^2.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*6.0+y.^2.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*6.0+y.^2.*(x+y-1.0).^2.*d2(:,3,:).*6.0+x.*y.*basis.*2.4e+1+x.^2.*y.^2.*d2(:,1,:).*6.0+y.*grad(:,1,:).*(x+y-1.0).^2.*1.2e+1+x.^2.*y.^2.*d2(:,3,:).*6.0+x.*y.^2.*grad(:,1,:).*2.4e+1+x.^2.*y.*grad(:,1,:).*1.2e+1+x.*y.^2.*grad(:,2,:).*1.2e+1+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d3(:,1,:)+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d3(:,3,:).*3.0+x.^2.*y.^2.*(x+y-1.0).^2.*d4(:,4,:)+x.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,1,:).*6.0+x.^2.*y.*(x.*2.0+y.*2.0-2.0).*d2(:,1,:).*6.0+x.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,3,:).*1.2e+1+x.^2.*y.*(x+y-1.0).^2.*d3(:,1,:).*2.0+x.*y.^2.*(x+y-1.0).^2.*d3(:,3,:).*6.0+x.*y.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*2.4e+1+x.*y.*(x+y-1.0).^2.*d2(:,1,:).*1.2e+1;
            val(:,5,:) = x.^2.*basis.*1.2e+1+x.*basis.*(x.*2.0+y.*2.0-2.0).*1.2e+1+x.^2.*grad(:,1,:).*(x.*2.0+y.*2.0-2.0).*6.0+x.^2.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*6.0+x.^2.*(x+y-1.0).^2.*d2(:,3,:).*6.0+x.*y.*basis.*2.4e+1+x.*grad(:,2,:).*(x+y-1.0).^2.*1.2e+1+x.^2.*y.^2.*d2(:,3,:).*6.0+x.^2.*y.^2.*d2(:,2,:).*6.0+x.^2.*y.*grad(:,1,:).*1.2e+1+x.*y.^2.*grad(:,2,:).*1.2e+1+x.^2.*y.*grad(:,2,:).*2.4e+1+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d3(:,4,:).*3.0+x.^2.*y.^2.*(x.*2.0+y.*2.0-2.0).*d3(:,2,:)+x.^2.*y.^2.*(x+y-1.0).^2.*d4(:,5,:)+x.^2.*y.*(x.*2.0+y.*2.0-2.0).*d2(:,3,:).*1.2e+1+x.*y.^2.*(x.*2.0+y.*2.0-2.0).*d2(:,2,:).*6.0+x.^2.*y.*(x.*2.0+y.*2.0-2.0).*d2(:,2,:).*6.0+x.^2.*y.*(x+y-1.0).^2.*d3(:,4,:).*6.0+x.*y.^2.*(x+y-1.0).^2.*d3(:,2,:).*2.0+x.*y.*grad(:,2,:).*(x.*2.0+y.*2.0-2.0).*2.4e+1+x.*y.*(x+y-1.0).^2.*d2(:,2,:).*1.2e+1;
        end
    end
end
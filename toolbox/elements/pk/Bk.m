classdef Bk < DGDk
    %BK - Cubic bubble times Dubiner polynomials of degree k-3
    %   OBJ = BK(DEGREE,MESH) creates the element bubble space of
    %   polynomial degree DEGREE >= 3 (a scalar or one value per element)
    %   on the triangular mesh MESH whose element basis functions are the
    %   Dubiner polynomials of DGDk of total degree DEGREE-3 multiplied by
    %   the cubic bubble b_T = x*y*(1-x-y) of the reference triangle.
    %   The basis functions vanish on the boundary of their element.
    %   Bk(3,MESH) is the cubic bubble of the MINI element.
    %
    %   BK functions:
    %       getElementBasis     - Evaluate the bubble element basis
    %       getGradElementBasis - Evaluate the gradient of the bubbles
    %       getD2ElementBasis   - Evaluate the Hessian of the bubbles
    %       getD3ElementBasis   - Evaluate third derivatives of bubbles
    %       getD4ElementBasis   - Evaluate fourth derivatives of bubbles
    %
    %   See also B2k, CurlBk, DGDk, Pk

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = Bk(degree,mesh)
            %BK - Create a cubic bubble element space
            if any(degree(:) < 3)
                error('fem4m:Bk:invalidDegree','Bk degree must be at least 3');
            end
            obj = obj@DGDk(degree-3,mesh);
            obj.degree4e = obj.degree4e+3;
            obj.degree = obj.degree+3;
            obj.degree4DoF = obj.degree4DoF+3;
        end
    end

    methods (Static)
        %% Elem basis helper functions
        function val = getElementBasis(p,x)
            %getElementBasis - Evaluate the bubble element basis
            %   VAL = getElementBasis(P,X) multiplies the Dubiner
            %   basis of total degree P-3 with the cubic bubble at
            %   the reference point X (1-by-2). VAL is
            %   1-by-(P-2)(P-1)/2.
            basis = DGDk.getElementBasis(p-3,x);
            y = x(:,2)';
            x = x(:,1)';
            val = (1-x-y).*x.*y.*basis;
        end

        function val = getGradElementBasis(p,x)
            %getGradElementBasis - Evaluate the gradient of the bubbles
            %   VAL = getGradElementBasis(P,X) evaluates the
            %   gradients of the cubic bubble times the Dubiner
            %   basis of total degree P-3 at the reference point X
            %   (1-by-2). VAL is (P-2)(P-1)/2-by-2.
            basis = DGDk.getElementBasis(p-3,x)';
            grad = DGDk.getGradElementBasis(p-3,x);
            y = x(:,2)';
            x = x(:,1)';
            val(:,1,:) = -y.*basis.*(x+y-1.0)-x.*y.*basis-x.*y.*grad(:,1,:).*(x+y-1.0);
            val(:,2,:) = -x.*basis.*(x+y-1.0)-x.*y.*basis-x.*y.*grad(:,2,:).*(x+y-1.0);
        end

        function val = getD2ElementBasis(p,x)
            %getD2ElementBasis - Evaluate the Hessian of the bubbles
            %   VAL = getD2ElementBasis(P,X) evaluates the Hessians
            %   [dxx,dyy,dxy] of the cubic bubble times the Dubiner
            %   basis of total degree P-3 at the reference point X
            %   (1-by-2). VAL is (P-2)(P-1)/2-by-3.
            basis = DGDk.getElementBasis(p-3,x)';
            grad = DGDk.getGradElementBasis(p-3,x);
            d2 = DGDk.getD2ElementBasis(p-3,x);
            y = x(:,2)';
            x = x(:,1)';
            val(:,1,:) = y.*basis.*-2.0-y.*grad(:,1,:).*(x+y-1.0).*2.0-x.*y.*grad(:,1,:).*2.0-x.*y.*(x+y-1.0).*d2(:,1,:);
            val(:,2,:) = x.*basis.*-2.0-x.*grad(:,2,:).*(x+y-1.0).*2.0-x.*y.*grad(:,2,:).*2.0-x.*y.*(x+y-1.0).*d2(:,2,:);
            val(:,3,:) = -basis.*(x+y-1.0)-x.*basis-y.*basis-x.*grad(:,1,:).*(x+y-1.0)-y.*grad(:,2,:).*(x+y-1.0)...
                -x.*y.*grad(:,1,:)-x.*y.*grad(:,2,:)-x.*y.*(x+y-1.0).*d2(:,3,:);
        end

        function val = getD3ElementBasis(p,x)
            %getD3ElementBasis - Evaluate third derivatives of bubbles
            %   VAL = getD3ElementBasis(P,X) evaluates the third
            %   derivatives [u_xxx,u_yyy,u_xxy,u_xyy] of the cubic
            %   bubble times the Dubiner basis of total degree P-3 at
            %   the reference point X (1-by-2). VAL is
            %   (P-2)(P-1)/2-by-4.
            basis = DGDk.getElementBasis(p-3,x)';
            grad = DGDk.getGradElementBasis(p-3,x);
            d2 = DGDk.getD2ElementBasis(p-3,x);
            d3 = DGDk.getD3ElementBasis(p-3,x);
            y = x(:,2)';
            x = x(:,1)';
            val(:,1,:) = y.*grad(:,1,:).*-6.0-y.*(x+y-1.0).*d2(:,1,:).*3.0-x.*y.*d2(:,1,:).*3.0-x.*y.*(x+y-1.0).*d3(:,1,:);
            val(:,2,:) = x.*grad(:,2,:).*-6.0-x.*(x+y-1.0).*d2(:,2,:).*3.0-x.*y.*d2(:,2,:).*3.0-x.*y.*(x+y-1.0).*d3(:,2,:);
            val(:,3,:) = x.*grad(:,1,:).*-2.0-y.*grad(:,1,:).*2.0-y.*grad(:,2,:).*2.0-basis.*2.0-grad(:,1,:).*(x+y-1.0).*2.0-x.*(x+y-1.0).*d2(:,1,:)-y.*(x+y-1.0).*d2(:,3,:).*2.0-x.*y.*d2(:,1,:)-x.*y.*d2(:,3,:).*2.0-x.*y.*(x+y-1.0).*d3(:,3,:);
            val(:,4,:) = x.*grad(:,1,:).*-2.0-x.*grad(:,2,:).*2.0-y.*grad(:,2,:).*2.0-basis.*2.0-grad(:,2,:).*(x+y-1.0).*2.0-x.*(x+y-1.0).*d2(:,3,:).*2.0-y.*(x+y-1.0).*d2(:,2,:)-x.*y.*d2(:,3,:).*2.0-x.*y.*d2(:,2,:)-x.*y.*(x+y-1.0).*d3(:,4,:);
        end

        function val = getD4ElementBasis(p,x)
            %getD4ElementBasis - Evaluate fourth derivatives of bubbles
            %   VAL = getD4ElementBasis(P,X) evaluates the fourth
            %   derivatives [u_xxxx,u_yyyy,u_xxyy,u_xxxy,u_xyyy] of
            %   the cubic bubble times the Dubiner basis of total
            %   degree P-3 at the reference point X (1-by-2). VAL is
            %   (P-2)(P-1)/2-by-5.
            grad = DGDk.getGradElementBasis(p-3,x);
            d2 = DGDk.getD2ElementBasis(p-3,x);
            d3 = DGDk.getD3ElementBasis(p-3,x);
            d4 = DGDk.getD4ElementBasis(p-3,x);
            y = x(:,2)';
            x = x(:,1)';
            val(:,1,:) = y.*d2(:,1,:).*-1.2e+1-y.*(x+y-1.0).*d3(:,1,:).*4.0-x.*y.*d3(:,1,:).*4.0-x.*y.*(x+y-1.0).*d4(:,1,:);
            val(:,2,:) = x.*d2(:,2,:).*-1.2e+1-x.*(x+y-1.0).*d3(:,2,:).*4.0-x.*y.*d3(:,2,:).*4.0-x.*y.*(x+y-1.0).*d4(:,2,:);
            val(:,3,:) = grad(:,1,:).*-4.0-grad(:,2,:).*4.0-(x+y-1.0).*d2(:,3,:).*4.0-x.*d2(:,1,:).*2.0-x.*d2(:,3,:).*4.0-y.*d2(:,3,:).*4.0-y.*d2(:,2,:).*2.0-x.*(x+y-1.0).*d3(:,3,:).*2.0-y.*(x+y-1.0).*d3(:,4,:).*2.0-x.*y.*d3(:,3,:).*2.0-x.*y.*d3(:,4,:).*2.0-x.*y.*(x+y-1.0).*d4(:,3,:);
            val(:,4,:) = grad(:,1,:).*-6.0-(x+y-1.0).*d2(:,1,:).*3.0-x.*d2(:,1,:).*3.0-y.*d2(:,1,:).*3.0-y.*d2(:,3,:).*6.0-x.*(x+y-1.0).*d3(:,1,:)-y.*(x+y-1.0).*d3(:,3,:).*3.0-x.*y.*d3(:,1,:)-x.*y.*d3(:,3,:).*3.0-x.*y.*(x+y-1.0).*d4(:,4,:);
            val(:,5,:) = grad(:,2,:).*-6.0-(x+y-1.0).*d2(:,2,:).*3.0-x.*d2(:,3,:).*6.0-x.*d2(:,2,:).*3.0-y.*d2(:,2,:).*3.0-x.*(x+y-1.0).*d3(:,4,:).*3.0-y.*(x+y-1.0).*d3(:,2,:)-x.*y.*d3(:,4,:).*3.0-x.*y.*d3(:,2,:)-x.*y.*(x+y-1.0).*d4(:,5,:);
        end
    end
end
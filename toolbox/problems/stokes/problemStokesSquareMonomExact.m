classdef problemStokesSquareMonomExact < ProblemStokesTemplateExact
    %problemStokesSquareMonomExact - Stokes problem with polynomial solution
    %   problemStokesSquareMonomExact is a Stokes problem
    %   -NU*Delta u + grad p = f, div u = 0 on the unit square with the
    %   divergence-free velocity u = curl(psi) = (d_y psi, -d_x psi) of the
    %   stream function psi = x^I y^J (a velocity of degree I+J-1), the
    %   zero-mean pressure p = x^M - 1/(M+1) and inhomogeneous Dirichlet
    %   data. It is meant for patch tests.
    %
    %   OBJ = problemStokesSquareMonomExact(NU,I,J,M) sets the viscosity NU,
    %   the stream function degrees I and J and the pressure degree M
    %   (M = 0 gives p = 0). All four arguments are required.
    %
    %   problemStokesSquareMonomExact functions:
    %       f           - Body force -NU*Delta u + grad p
    %       exactU      - Exact velocity
    %       exactGradU  - Velocity gradient
    %       exactP      - Exact pressure
    %
    %   problemStokesSquareMonomExact properties:
    %       geom        - Unit square geometry
    %       i           - Stream function degree in x
    %       j           - Stream function degree in y
    %       m           - Pressure degree
    %
    %   See also ProblemStokesTemplateExact, problemStokesExactSquareP3Sigma

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        i
        j
        m
    end

    methods

        function obj = problemStokesSquareMonomExact(nu,i,j,m)
            obj = obj@ProblemStokesTemplateExact(nu);
            obj.i = i;
            obj.j = j;
            obj.m = m;
        end

        function val = f(obj,x)
            val = -obj.nu*[obj.psi(x,2,1) + obj.psi(x,0,3),...
                -obj.psi(x,3,0) - obj.psi(x,1,2)];
            val(:,1) = val(:,1) + monomial(x,obj.m,0,1,0);
        end

        function val = exactU(obj,x)
            val = [obj.psi(x,0,1),-obj.psi(x,1,0)];
        end

        function val = exactGradU(obj,x)
            val = zeros(size(x,1),2,2);
            val(:,1,1) = obj.psi(x,1,1);
            val(:,1,2) = obj.psi(x,0,2);
            val(:,2,1) = -obj.psi(x,2,0);
            val(:,2,2) = -obj.psi(x,1,1);
        end

        function val = exactP(obj,x)
            if obj.m == 0
                val = zeros(size(x,1),1);
            else
                val = monomial(x,obj.m,0,0,0) - 1/(obj.m+1);
            end
        end

    end

    methods (Access = private)

        function val = psi(obj,x,dx,dy)
            val = monomial(x,obj.i,obj.j,dx,dy);
        end

    end

end

function val = monomial(x,i,j,dx,dy)
%MONOMIAL - Derivative d^(dx,dy) of x^i y^j
val = derivative1d(x(:,1),i,dx).*derivative1d(x(:,2),j,dy);
end

function val = derivative1d(x,k,d)
%derivative1d - d-th derivative of x^k
if d > k
    val = zeros(size(x));
else
    val = prod(k-d+1:k)*x.^(k-d);
end
end
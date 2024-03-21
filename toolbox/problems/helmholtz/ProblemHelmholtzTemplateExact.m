classdef ProblemHelmholtzTemplateExact
    %ProblemHelmholtzTemplateExact - Abstract Helmholtz problem with an exact solution
    %   ProblemHelmholtzTemplateExact is an abstract template for Helmholtz
    %   problems -Laplace(u) - omega^2*refraction(x)^2*u = f with a known
    %   exact solution. Subclasses supply the geometry geom and the exact
    %   solution (exactU, exactGradU); the template derives the Dirichlet
    %   data u = u4Db, the Neumann data grad(u)*n = g4Nb and the impedance
    %   data grad(u)*n + 1i*omega*refraction(x)*u = g4Rb from it. Solvers
    %   impose g4Rb on the Neumann boundary (dsn terms) and u4Db on the
    %   Dirichlet boundary. The source term is f = 0, so
    %   the exact solution must solve the homogeneous equation unless the
    %   subclass overrides f.
    %
    %   OBJ = ProblemHelmholtzTemplateExact(OMEGA) sets the angular
    %   frequency OMEGA (the wavenumber where the refractive index is 1);
    %   the subclass constructors call it. A missing or empty OMEGA raises
    %   the error fem4m:ProblemHelmholtzTemplateExact:missingOmega.
    %
    %   The data functions are evaluated at the points X (N-by-d, one row
    %   per point) and return an N-by-1 vector; exactGradU returns the
    %   N-by-d gradient (one row per point), and the normals passed to g4Nb
    %   and g4Rb are the outward unit normals (N-by-d, one row per point).
    %
    %   ProblemHelmholtzTemplateExact functions:
    %       f           - Zero source term
    %       u4Db        - Dirichlet data, taken from exactU
    %       g4Nb        - Neumann data grad(u)*n, from exactGradU
    %       g4Rb        - Impedance data grad(u)*n + 1i*omega*refraction*u
    %       refraction  - Refractive index (default: 1)
    %       exactU      - Exact solution (abstract)
    %       exactGradU  - Exact gradient of the solution (abstract)
    %
    %   ProblemHelmholtzTemplateExact properties:
    %       geom        - Domain geometry (abstract)
    %       omega       - Angular frequency (wavenumber for refraction 1)
    %
    %   See also problemHelmholtzSquareExact, problemHelmholtzSquareExact2,
    %   problemHelmholtzPlaneWave, problemHelmholtzLshapeExact,
    %   problemHelmholtzLshapeSmoothExact, problemHelmholtzSquareReflection,
    %   problemHelmholtzGaussianBeam

    % Copyright (c) 2026 by the fem4m authors

    properties(Abstract)
        geom
    end

    properties
        omega
    end

    methods(Abstract)
        val = exactU(x);
        val = exactGradU(x);
    end

    methods

        function obj = ProblemHelmholtzTemplateExact(omega)
            if nargin < 1 || isempty(omega)
                error('fem4m:ProblemHelmholtzTemplateExact:missingOmega','%s: angular frequency omega required',class(obj));
            end
            obj.omega = omega;
        end

        function val = f(~,x)
            val = zeros(size(x,1),1);
        end

        function val = u4Db(obj,x)
            val = obj.exactU(x);
        end

        function val = g4Nb(obj,x,normals4Nb)
            val = sum(obj.exactGradU(x).*normals4Nb,2);
        end

        function val = g4Rb(obj,x,normals4Sb)
            val = sum(obj.exactGradU(x).*normals4Sb,2) + 1i*obj.omega.*obj.refraction(x).*obj.exactU(x);
        end

        function val = refraction(~,x)
            val = ones(size(x,1),1);
        end

    end

end
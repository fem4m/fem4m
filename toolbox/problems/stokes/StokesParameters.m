classdef StokesParameters
    %StokesParameters - Material parameters for Stokes flow
    %   OBJ = StokesParameters() creates Stokes parameters with the
    %   viscosity NU = 1.
    %
    %   OBJ = StokesParameters(NU) sets the viscosity NU; an empty NU gives
    %   NU = 1.
    %
    %   The matrices A, AWeakSym and A3d are computed from NU on access.
    %   They map a stress sigma to dev(sigma)/(2*NU); A and AWeakSym are
    %   the limits of ElasticityParameters.CInv and ElasticityParameters.A
    %   for lambda -> inf and mu = NU.
    %
    %   StokesParameters properties:
    %       nu          - Viscosity
    %       A           - Compliance matrix on symmetric 2d tensors stored
    %                     as [11 22 12]
    %       AWeakSym    - Compliance matrix on general 2d tensors stored
    %                     column-wise as [11 21 12 22]
    %       A3d         - Compliance matrix on general 3d tensors (9
    %                     entries, the diagonal at 1, 5 and 9)
    %
    %   See also ProblemStokesTemplate, ElasticityParameters

    % Copyright (c) 2026 by the fem4m authors

    properties
        nu
        A
        AWeakSym
        A3d
    end

    methods

        function obj = StokesParameters(nu)
            if nargin < 1 || isempty(nu)
                nu = 1;
            end
            obj.nu = nu;
        end

        function A = get.A(obj)
            A  = 1/(4*obj.nu)*[1 -1 0; -1 1 0 ; 0 0 2];
        end

        function A = get.AWeakSym(obj)
            A  = 1/(4*obj.nu)*[1 0 0 -1 ; 0 2 0 0 ; 0 0 2 0; -1 0 0 1 ];
        end

        function A = get.A3d(obj)
            A = eye(9);
            A(1:4:9,1:4:9) = A(1:4:9,1:4:9)-1/3;
            A  = 1/(2*obj.nu)*A;
        end

    end

end
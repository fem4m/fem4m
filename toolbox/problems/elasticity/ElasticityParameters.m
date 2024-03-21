classdef ElasticityParameters
    %ElasticityParameters - Material parameters for linear elasticity
    %   OBJ = ElasticityParameters(NAME1,VALUE1,NAME2,VALUE2) creates
    %   material constants from the pair (E, nu) or the pair (lambda, mu)
    %   specified by NAME1 and NAME2, in either order and ignoring case.
    %   The other pair is derived: lambda = E*nu/((1+nu)*(1-2*nu)) and
    %   mu = E/(2*(1+nu)), or E = mu*(3*lambda+2*mu)/(lambda+mu) and
    %   nu = lambda/(2*(lambda+mu)). Fewer than four arguments or an
    %   empty NAME1 or NAME2 raise the error
    %   fem4m:ElasticityParameters:missingParameters, any other pair of
    %   names fem4m:ElasticityParameters:unknownParameterPair.
    %
    %   The 2d tensors C and CInv act on symmetric tensors stored as
    %   [11 22 12] (strain [e11 e22 e12], stress [s11 s22 s12], without
    %   the factor 2 of the engineering shear strain); A and AInv act on
    %   general 2-by-2 tensors stored column-wise as [11 21 12 22].
    %
    %   ElasticityParameters functions:
    %       M           - Bending stiffness matrix of a plate
    %       D           - Flexural rigidity of a plate
    %
    %   ElasticityParameters properties:
    %       lambda      - First Lamé parameter
    %       mu          - Shear modulus (second Lamé parameter)
    %       E           - Young's modulus
    %       nu          - Poisson ratio
    %       C           - Stiffness matrix, maps [e11 e22 e12] to
    %                     [s11 s22 s12]
    %       CInv        - Compliance matrix, inverse of C
    %       A           - Compliance tensor on the components
    %                     [11 21 12 22] of a 2-by-2 stress
    %       AInv        - Stiffness tensor on the components
    %                     [11 21 12 22], inverse of A
    %
    %   See also ProblemElasticityTemplate, ProblemElasticityTemplateExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        lambda
        mu
        E
        nu
        C
        CInv
        A
        AInv
    end

    methods

        function obj = ElasticityParameters(name1,value1,name2,value2)
            if nargin < 4 || isempty(name1) || isempty(name2)
                error('fem4m:ElasticityParameters:missingParameters','%s: material parameters required, call with (''E'',E,''nu'',nu) or (''lambda'',lambda,''mu'',mu)',class(obj));
            end
            if strcmpi(name1,'E'), E = value1; elseif strcmpi(name2,'E'), E = value2; else, E = []; end
            if strcmpi(name1,'nu'), nu = value1; elseif strcmpi(name2,'nu'), nu = value2; else, nu = []; end
            if strcmpi(name1,'lambda'), lambda = value1; elseif strcmpi(name2,'lambda'), lambda = value2; else, lambda = []; end
            if strcmpi(name1,'mu'), mu = value1; elseif strcmpi(name2,'mu'), mu = value2; else, mu = []; end
            if ~isempty(E) && ~isempty(nu)
                obj.lambda = (E*nu)/((1+nu)*(1-2*nu));
                obj.mu = E/(2*(1+nu));
                obj.E = E;
                obj.nu = nu;
            elseif ~isempty(lambda) && ~isempty(mu)
                obj.lambda = lambda;
                obj.mu = mu;
                obj.E = mu*(3*lambda+2*mu)/(lambda+mu);
                obj.nu = lambda/(2*(lambda+mu));
            else
                error('fem4m:ElasticityParameters:unknownParameterPair','Elasticity parameter pair (%s,%s) not implemented',name1,name2)
            end
        end

        function C = get.C(obj)
            C = [2*obj.mu+obj.lambda,obj.lambda,0;obj.lambda,2*obj.mu+obj.lambda,0;0,0,2*obj.mu];
        end

        function CInv = get.CInv(obj)
            CInv = (obj.lambda+2*obj.mu)/(4*obj.mu*(obj.lambda+obj.mu))*[1 0 0; 0 1 0 ; 0 0 0] +...
                -obj.lambda/(4*obj.mu*(obj.lambda+obj.mu))*[ 0 1 0; 1 0 0; 0 0 0] +...
                1/(2*obj.mu)*[0 0 0 ; 0 0 0 ; 0 0 1];
        end

        function AInv = get.AInv(obj)
            AInv = [2*obj.mu+obj.lambda,0,0,obj.lambda;0,2*obj.mu,0,0;0,0,2*obj.mu,0;obj.lambda,0,0,2*obj.mu+obj.lambda];
        end

        function A = get.A(obj)
            A  = (obj.lambda+2*obj.mu)/(4*obj.mu*(obj.lambda+obj.mu))*[1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1 ] +...
                -obj.lambda/(4*obj.mu*(obj.lambda+obj.mu))*[ 0 0 0 1 ; 0 0 0 0; 0 0 0 0; 1 0 0 0 ] +...
                1/(2*obj.mu)*[0 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 0 ];
        end

        function M = M(obj,t)
            %M - Bending stiffness matrix of a plate
            %   VAL = M(OBJ,T) returns D(OBJ,T)*[1 nu 0;nu 1 0;0 0 1-nu]
            %   for the plate thickness T.
            M = obj.D(t)*[1,obj.nu,0; obj.nu,1,0; 0,0, 1-obj.nu];
        end

        function D = D(obj,t)
            %D - Flexural rigidity of a plate
            %   VAL = D(OBJ,T) returns E*T^3/(12*(1-nu^2)) for the plate
            %   thickness T.
            D = t^3*obj.E/12/(1-obj.nu^2);
        end
    end

end
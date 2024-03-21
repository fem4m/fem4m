classdef problemElasticityCrack < ProblemElasticityTemplate
    %problemElasticityCrack - Elasticity on a cracked (slit) domain
    %   OBJ = problemElasticityCrack(NAME1,VALUE1,NAME2,VALUE2) supplies
    %   a linear elasticity problem on the slit domain [-1,1]^2 (slit from
    %   (0,0) to (1,0)) with zero body force, zero displacement data and
    %   the vertical traction (0,1) on y = 1 and (0,-1) on y = -1. The
    %   traction acts only on Neumann parts of geom; the default geometry
    %   is Dirichlet on the whole boundary, so the solution is zero. The
    %   material parameters NAME1,VALUE1,NAME2,VALUE2 are given as in
    %   ElasticityParameters and are required.
    %
    %   problemElasticityCrack functions:
    %       f     - Zero body force
    %       u4Db  - Zero displacement data
    %       g     - Vertical traction +1 on y = 1 and -1 on y = -1
    %
    %   problemElasticityCrack properties:
    %       geom  - Slit domain geometry (Dirichlet boundary)
    %       scale - Amplification of the displacement plots
    %
    %   See also ProblemElasticityTemplate, ElasticityParameters,
    %   problemElasticityCooks

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Slit'))
        scale = 1000
    end


    methods

        function val = f(~,x)
            val = zeros(size(x,1),2);
        end

        function val = u4Db(~,x)
            val = zeros(size(x,1),2);
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),2);
            val(x(:,2) == 1,2) = 1;
            val(x(:,2) == -1,2) = -1;
        end

    end
end
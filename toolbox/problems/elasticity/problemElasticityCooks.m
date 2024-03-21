classdef problemElasticityCooks < ProblemElasticityTemplate
    %problemElasticityCooks - Linear elasticity of the Cook membrane
    %   OBJ = problemElasticityCooks(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies a linear elasticity problem on the Cook membrane, the
    %   quadrilateral with the corners (0,-44), (48,0), (48,16) and (0,0),
    %   with zero body force, zero displacement on the clamped left edge
    %   x = 0, a unit vertical traction on the right edge x = 48 and zero
    %   traction on the other edges. The material parameters
    %   NAME1,VALUE1,NAME2,VALUE2 are given as in ElasticityParameters and
    %   are required.
    %
    %   problemElasticityCooks functions:
    %       f          - Zero body force
    %       u4Db       - Zero displacement data
    %       gradu4Dbt  - Tangential derivative grad(u)*t (zero)
    %       gradu4Dbtt - Tangential strain t'*grad(u)*t (zero)
    %       d2u4Dbntt  - Second tangential derivative n'*(d_tt u) (zero)
    %       g          - Traction sigma*n: (0,1) on x = 48, zero elsewhere
    %
    %   problemElasticityCooks properties:
    %       geom       - Cook membrane geometry (Cooks_*.dat)
    %       scale      - Amplification of the displacement plots
    %
    %   See also ProblemElasticityTemplate, ElasticityParameters,
    %   problemElasticityCrack, problemElasticityEigenvalueCooks

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cooks'))
        scale = 2000
    end


    methods

        function val = f(~,x)
            val = zeros(size(x,1),2);
        end

        function val = u4Db(~,x)
            val = zeros(size(x,1),2);
        end

        function val = gradu4Dbt(~,x,~)
            val = zeros(size(x,1),2);
        end

        function val = gradu4Dbtt(~,x,~)
            val = zeros(size(x,1),1);
        end

        function val = d2u4Dbntt(~,x,~)
            val = zeros(size(x,1),1);
        end

        function val = g(~,x,normals)
            val = zeros(size(x,1),2);
            val(x(:,1) == 48 & normals(:,2)==0,2) = 1;
        end

    end
end
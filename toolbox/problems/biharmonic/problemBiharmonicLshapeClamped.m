classdef problemBiharmonicLshapeClamped < ProblemBiharmonicTemplate
    %problemBiharmonicLshapeClamped - Clamped biharmonic on an L-shape
    %   problemBiharmonicLshapeClamped is a concrete clamped biharmonic
    %   problem on the L-shape (-1,1)^2 \ [0,1]x[-1,0] with the constant
    %   source term f = 1 and homogeneous boundary data.
    %
    %   problemBiharmonicLshapeClamped functions:
    %       f     - Constant source term 1
    %       g     - Zero normal-derivative data
    %       u4Db  - Zero Dirichlet data
    %
    %   problemBiharmonicLshapeClamped properties:
    %       geom  - L-shaped domain geometry
    %
    %   See also ProblemBiharmonicTemplate, problemBiharmonicLshapeClampedExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape'))
    end


    methods

        function val = f(~,x,varargin)
            val = ones(size(x,1),1);
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),1);
        end

        function val = u4Db(~,x)
            val = zeros(size(x,1),1);
        end

    end

end
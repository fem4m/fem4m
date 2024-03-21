classdef problemBendingSquare < ProblemBiharmonicTemplate
    %problemBendingSquare - Nonlinear plate bending on a unit square
    %   problemBendingSquare is a bending problem for nonlinear plate
    %   models on the unit square, clamped (Dirichlet) on the left and
    %   right edges and free (Neumann) on the bottom and top edges, with
    %   a constant transverse load and zero boundary data. All data have
    %   three components (N-by-3); the load 1e-2 acts on the third
    %   component.
    %
    %   problemBendingSquare functions:
    %       f     - Constant transverse load (third component)
    %       u4Db  - Zero Dirichlet data
    %       g     - Zero normal-derivative data
    %
    %   problemBendingSquare properties:
    %       geom  - Square geometry
    %
    %   See also ProblemBiharmonicTemplate, problemBiharmonicLshapeClamped

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry',['N','D','N','D'],'reftype','cross'))
    end


    methods

        function val = f(~,x)
            val = zeros(size(x,1),3);
            val(:,3) = 1e-2;
        end

        function val = u4Db(~,x)
            val = zeros(size(x,1),3);
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),3);
        end

    end

end
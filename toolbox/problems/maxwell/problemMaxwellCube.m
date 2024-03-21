classdef problemMaxwellCube
    %problemMaxwellCube - Maxwell on a cube, piecewise-constant permeability
    %   problemMaxwellCube is a Maxwell source problem
    %   curl(mu^(-1) curl(u)) = f on the unit cube with the constant source
    %   current f = (1,0,0), a piecewise-constant permeability mu and zero
    %   Dirichlet data. No exact solution is known.
    %
    %   OBJ = problemMaxwellCube() uses unit permeability.
    %
    %   OBJ = problemMaxwellCube(MU) sets the permeability to MU on the
    %   L-shaped prism outside the edge region y < 1/2, z < 1/2 along the
    %   x-axis and to 1 inside that region.
    %
    %   problemMaxwellCube functions:
    %       f     - Constant source current in x-direction
    %       curlf - Curl of the source current (zero)
    %       u4Db  - Zero Dirichlet field
    %       mu    - Piecewise-constant permeability
    %
    %   problemMaxwellCube properties:
    %       geom  - Cube geometry
    %       muVal - Permeability outside the edge region (default 1)
    %
    %   See also problemMaxwellCubeExact, problemMaxwellCubeExact2,
    %   problemMaxwellLshape3dExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1))
        muVal = 1;
    end


    methods

        function obj = problemMaxwellCube(mu)
            if nargin < 1 || isempty(mu)
                obj.muVal = 1;
            else
                obj.muVal = mu;
            end
        end

        function val = f(~,x)
            val = zeros(size(x,1),3);
            %val = ones(size(x,1),3);
            val(:,1) = 1;

        end

        function val = curlf(~,x)
            val = zeros(size(x,1),3);
        end

        function val = u4Db(~,x)
            val = zeros(size(x,1),3);
        end

        function val = mu(obj,x)
            %MU - Piecewise-constant permeability
            val = ones(size(x,1),1);
            % checkerboard 3d
            %val(x(:,1)>0.5 & x(:,2)>0.5 & x(:,3)>0.5) = obj.muVal;
            %val(x(:,1)<0.5 & x(:,2)<0.5 & x(:,3)>0.5) = obj.muVal;
            %val(x(:,1)>0.5 & x(:,2)<0.5 & x(:,3)<0.5) = obj.muVal;
            %val(x(:,1)<0.5 & x(:,2)>0.5 & x(:,3)<0.5) = obj.muVal;
            % checkerboard 2d perp
            %val(x(:,1)<0.5 & x(:,2)<0.5) = obj.muVal;
            %val(x(:,1)>0.5 & x(:,2)>0.5) = obj.muVal;
            % checkerboard 2d parallel
            %val(x(:,3)<0.5 & x(:,2)<0.5) = obj.muVal;
            %val(x(:,3)>0.5 & x(:,2)>0.5) = obj.muVal;
            % Lshaped parallel
            val(~(x(:,3)<0.5 & x(:,2)<0.5)) = obj.muVal;
        end


    end
end
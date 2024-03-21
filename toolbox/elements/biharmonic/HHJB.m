classdef HHJB < HHJk
    %HHJB - Hellan-Herrmann-Johnson element of degree k with bubbles
    %   OBJ = HHJB(DEGREE,MESH) creates the Hellan-Herrmann-Johnson
    %   element HHJk(DEGREE,MESH) enriched by 2*DEGREE divdiv-free
    %   element bubbles of degree DEGREE+1 per element. The TDNNS
    %   (tangential-displacement normal-normal-stress) method for linear
    %   elasticity uses this space; plate problems do not need the
    %   bubbles. HHJB(0,MESH) has no bubbles and is the space of
    %   HHJk(0,MESH). HHJB has the functions and properties of HHJk,
    %   with bubbles = true.
    %
    %   See also HHJk, FE

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = HHJB(degree,mesh)
            %HHJB - Create a Hellan-Herrmann-Johnson element with bubbles
            obj = obj@HHJk(degree,mesh);
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) switches the bubbles on and sets the
            %   numbers of HHJk.dofCounts, which then include the 2*DEGREE
            %   interior bubble dofs.
            obj.bubbles = true;
            obj = dofCounts@HHJk(obj);
        end

    end
end
classdef HBDMk < BDMk
    %HBDMk - Hybridized Brezzi-Douglas-Marini element of degree k
    %   OBJ = HBDMk(DEGREE,MESH) creates the hybridized variant of the
    %   Brezzi-Douglas-Marini element of polynomial degree DEGREE >= 1
    %   on the triangular mesh MESH, where all degrees of freedom are
    %   stored element-locally and continuity is imposed externally by
    %   Lagrange multipliers. HBDMk has the functions and properties of
    %   BDMk; as it has no edge dofs, dof4Nb and fixedNodes are empty
    %   and projectNbBoundary returns zeros.
    %
    %   See also BDMk, HRTk, RTk

    % Copyright (c) 2026 by the fem4m authors

    methods (Access = protected)
        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) treats all edge and element degrees of
            %   freedom of the hybridized Brezzi-Douglas-Marini element as
            %   element-local ones.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 0;
            obj.nrBasis4e = obj.nrBasisFunc;
        end
    end

end
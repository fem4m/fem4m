classdef MarkP < Mark
    %MarkP - Pure p-refinement marking strategy
    %   OBJ = MarkP() creates a pure p-refinement marking strategy. It
    %   has no parameters; the arguments of the Mark constructor are
    %   accepted and not used.
    %
    %   MarkP functions:
    %       mark - Mark all elements for p-refinement
    %
    %   See also Mark, MarkHP, MarkFixedFraction

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = mark(obj,fe,varargin)
            %MARK - Mark all elements for p-refinement
            %   OBJ = MARK(OBJ,FE) marks every element of FE.mesh for
            %   p-refinement in OBJ.markedP4e and no edge or element for
            %   h-refinement.
            %
            %   OBJ = MARK(OBJ,FE,INDICATORS,...) ignores the indicator
            %   vectors.
            mesh = fe.mesh;
            obj = obj.initializeMark(mesh);
            obj.markedP4e = ones(mesh.nrElems,1);
        end
    end
end
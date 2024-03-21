classdef MarkUniform < Mark
    %MarkUniform - Uniform marking strategy
    %   OBJ = MarkUniform() creates a uniform marking strategy. It has
    %   no parameters; the arguments of the Mark constructor are
    %   accepted and not used.
    %
    %   MarkUniform functions:
    %       mark - Mark all edges and elements uniformly
    %
    %   See also Mark, MarkBulk, MarkMaximum

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = mark(obj,fe,varargin)
            %MARK - Mark all edges and elements uniformly
            %   OBJ = MARK(OBJ,FE) marks every edge and element of FE.mesh
            %   for h-refinement in OBJ.marked4ed and OBJ.marked4e and no
            %   element for p-refinement.
            %
            %   OBJ = MARK(OBJ,FE,INDICATORS,...) ignores the indicator
            %   vectors.
            mesh = fe.mesh;
            obj = obj.initializeMark(mesh);
            obj.marked4ed = ones(mesh.nrEdges,1);
            obj.marked4e = ones(mesh.nrElems,1);
        end
    end
end
classdef MarkMean < Mark
    %MarkMean - Mean marking strategy of the indicators
    %   OBJ = MarkMean() creates a mean marking strategy. It has no
    %   parameters; the arguments of the Mark constructor are accepted
    %   and not used.
    %
    %   MarkMean functions:
    %       mark - Mark edges and elements by the mean criterion
    %
    %   See also Mark, MarkBulk, MarkMaximum

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = mark(obj,fe,varargin)
            %MARK - Mark edges and elements by the mean criterion
            %   OBJ = MARK(OBJ,FE,INDICATORS) marks every entry of
            %   INDICATORS whose square is at least three quarters of the
            %   sum of all squared entries divided by the number of
            %   elements of FE.mesh, i.e. of the mean squared indicator
            %   for indicators per element, and stores the marks in
            %   OBJ.marked4ed and OBJ.marked4e.
            %   INDICATORS is a column vector given per edge, per element
            %   or per node. An element indicator also marks the edges of
            %   the element (not on polygonal meshes), a node indicator
            %   all edges and elements at the node. The kind is
            %   recognised by the length, checked in the order edges,
            %   elements, nodes; a vector of another length marks
            %   nothing.
            %
            %   OBJ = MARK(OBJ,FE,INDICATORS1,INDICATORS2,...) sums the
            %   squares of the entries of all indicator vectors.
            mesh = fe.mesh;
            obj = obj.initializeMark(mesh);
            N = length(varargin);
            indicators = [];
            dim = zeros(N+1,1);
            for k = 1 : N
                indicators = [indicators;varargin{k}]; %#ok<AGROW>
                dim(k+1) = length(indicators);
            end
            I = find(indicators.^2 >= 0.75*sum(indicators.^2)/mesh.nrElems);
            for k = 1 : N
                Ik = I(dim(k)< I & I<=dim(k+1));
                if length(varargin{k}) == size(mesh.n4ed,1)
                    obj.marked4ed(Ik-dim(k)) = 1;
                elseif length(varargin{k}) == size(mesh.n4e,1)
                    if ~iscell(mesh.ed4e)
                        obj.marked4ed(mesh.ed4e(Ik-dim(k),:)) = 1;
                    end
                    obj.marked4e(Ik-dim(k)) = 1;
                elseif length(varargin{k}) == size(mesh.c4n,1)
                    obj.marked4ed(nonzeros(mesh.ed4n(Ik-dim(k),:))) = 1;
                    obj.marked4e(nonzeros(mesh.e4n(Ik-dim(k),:))) = 1;
                end
            end
        end
    end
end
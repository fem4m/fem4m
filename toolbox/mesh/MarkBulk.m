classdef MarkBulk < Mark
    %MarkBulk - Bulk marking strategy with parameter theta
    %   OBJ = MarkBulk() creates a bulk marking strategy with the bulk
    %   parameter THETA = 0.5.
    %
    %   OBJ = MarkBulk(THETA) sets the bulk parameter THETA. The further
    %   arguments of the Mark constructor are accepted and not used.
    %
    %   MarkBulk functions:
    %       mark - Mark edges and elements by the bulk criterion
    %
    %   See also Mark, MarkMaximum, MarkMean

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = mark(obj,fe,varargin)
            %MARK - Mark edges and elements by the bulk criterion
            %   OBJ = MARK(OBJ,FE,INDICATORS) marks the largest entries of
            %   INDICATORS, as few as possible, whose squares sum to at
            %   least OBJ.theta times the sum of all squared entries, and
            %   stores the marks in OBJ.marked4ed and OBJ.marked4e of the
            %   mesh FE.mesh. For OBJ.theta >= 1 every entry is marked.
            %   INDICATORS is a column vector given per edge, per element
            %   or per node. An element indicator also marks the edges of
            %   the element (not on polygonal meshes), a node indicator
            %   all edges and elements at the node. The kind is
            %   recognised by the length, checked in the order edges,
            %   elements, nodes; a vector of another length marks
            %   nothing.
            %
            %   OBJ = MARK(OBJ,FE,INDICATORS1,INDICATORS2,...) applies the
            %   criterion to all entries of the indicator vectors
            %   together.
            mesh = fe.mesh;
            obj = obj.initializeMark(mesh);
            N = length(varargin);
            indicators = [];
            dim = zeros(N+1,1);
            for k = 1 : N
                indicators = [indicators;varargin{k}]; %#ok<AGROW>
                dim(k+1) = length(indicators);
            end
            [indicators,ind] = sort(indicators,'descend');
            J = find(cumsum(indicators.^2) >= obj.theta*sum(indicators.^2),1,'first');
            if obj.theta >= 1
                I = ind;
            else
                I = ind(1:J);
            end
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
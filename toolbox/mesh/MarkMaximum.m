classdef MarkMaximum < Mark
    %MarkMaximum - Maximum marking strategy with parameter theta
    %   OBJ = MarkMaximum() creates a maximum marking strategy with the
    %   parameter THETA = 0.5.
    %
    %   OBJ = MarkMaximum(THETA) sets the parameter THETA. The further
    %   arguments of the Mark constructor are accepted and not used.
    %
    %   MarkMaximum functions:
    %       mark - Mark edges and elements by the maximum criterion
    %
    %   See also Mark, MarkBulk, MarkMean

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = mark(obj,fe,varargin)
            %MARK - Mark edges and elements by the maximum criterion
            %   OBJ = MARK(OBJ,FE,INDICATORS) marks every entry of
            %   INDICATORS greater than OBJ.theta times the largest entry
            %   and stores the marks in OBJ.marked4ed and OBJ.marked4e of
            %   the mesh FE.mesh. INDICATORS is a column vector given per
            %   edge, per element or per node. An element indicator also
            %   marks the edges of the element (not on polygonal meshes),
            %   a node indicator all edges and elements at the node. The
            %   kind is recognised by the length, checked in the order
            %   edges, elements, nodes; a vector of another length marks
            %   nothing.
            %
            %   OBJ = MARK(OBJ,FE,INDICATORS1,INDICATORS2,...) compares
            %   all entries of the indicator vectors with the largest
            %   entry of all of them.
            mesh = fe.mesh;
            obj = obj.initializeMark(mesh);
            N = length(varargin);
            indicators = [];
            dim = zeros(N+1,1);
            for k = 1 : N
                indicators = [indicators;varargin{k}]; %#ok<AGROW>
                dim(k+1) = length(indicators);
            end
            I = find(indicators > obj.theta*max(indicators));
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
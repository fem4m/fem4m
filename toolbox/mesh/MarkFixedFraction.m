classdef MarkFixedFraction < Mark
    %MarkFixedFraction - Fixed fraction marking strategy scaled by the polynomial degree
    %   OBJ = MarkFixedFraction() creates a fixed fraction marking
    %   strategy. It has no parameters; the arguments of the Mark
    %   constructor are accepted and not used.
    %
    %   MarkFixedFraction functions:
    %       mark - Mark a fixed number of the largest indicators
    %
    %   See also Mark, MarkHP, MarkP

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = mark(obj,fe,varargin)
            %MARK - Mark a fixed number of the largest indicators
            %   OBJ = MARK(OBJ,FE,INDICATORS) marks the
            %   ceil(nrElems/(4*FE.degree)) largest entries of
            %   INDICATORS, nrElems the number of elements of FE.mesh and
            %   FE.degree the polynomial degree of FE (the maximum degree
            %   for HPk), and stores the marks in OBJ.marked4ed and
            %   OBJ.marked4e. INDICATORS must have at least that many
            %   entries.
            %   INDICATORS is a column vector given per edge, per element
            %   or per node. An element indicator also marks the edges of
            %   the element (not on polygonal meshes), a node indicator
            %   all edges and elements at the node. The kind is
            %   recognised by the length, checked in the order edges,
            %   elements, nodes; a vector of another length marks
            %   nothing.
            %
            %   OBJ = MARK(OBJ,FE,INDICATORS1,INDICATORS2,...) marks the
            %   largest entries of all indicator vectors together.
            mesh = fe.mesh;
            obj = obj.initializeMark(mesh);
            N = length(varargin);
            indicators = [];
            dim = zeros(N+1,1);
            for k = 1 : N
                indicators = [indicators;varargin{k}]; %#ok<AGROW>
                dim(k+1) = length(indicators);
            end
            [~,I] = sort(indicators,'descend');
            I = I(1:ceil(0.25/fe.degree*mesh.nrElems));
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
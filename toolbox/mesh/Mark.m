classdef Mark
    %MARK - Base class of mesh marking strategies for adaptive refinement
    %   Mark is abstract. Its subclasses mark edges and elements for
    %   h-refinement (MarkBulk, MarkMaximum, MarkMean, MarkFixedFraction,
    %   MarkUniform), elements for p-refinement (MarkP), or split the
    %   marked elements into h- and p-refinement (MarkHP, MarkHP3Edges).
    %   They inherit the constructor described below; an omitted or
    %   empty argument keeps its default. THETA is read by MarkBulk and
    %   MarkMaximum, GAMMAH, GAMMAP, GAMMAN and INITIALPREDICT by MarkHP
    %   and MarkHP3Edges; the other subclasses ignore all arguments.
    %
    %   OBJ = MarkBulk() creates a marking strategy with the default
    %   parameters THETA = 0.5, GAMMAH = 4, GAMMAP = 0.4, GAMMAN = 1 and
    %   INITIALPREDICT = inf.
    %
    %   OBJ = MarkBulk(THETA) also sets the bulk parameter THETA.
    %
    %   OBJ = MarkHP(THETA,GAMMAH) also sets the update parameter
    %   GAMMAH of the predictor for h-refined elements; a larger GAMMAH
    %   makes p-refinement of their children more likely. MarkHP and
    %   MarkHP3Edges ignore THETA.
    %
    %   OBJ = MarkHP(...,GAMMAP) also sets the update parameter GAMMAP
    %   of the predictor for p-refined elements; a smaller GAMMAP makes
    %   their next refinement more likely an h-refinement.
    %
    %   OBJ = MarkHP(...,GAMMAN) also sets the update parameter GAMMAN
    %   of the predictor for elements that are not refined; a smaller
    %   GAMMAN makes their next refinement more likely an h-refinement.
    %
    %   OBJ = MarkHP(...,INITIALPREDICT) also sets the initial value
    %   INITIALPREDICT of the predictor predict4e; the default inf makes
    %   the hp strategies p-refine every marked element first.
    %
    %   A subclass method mark(OBJ,FE,...) returns OBJ with the marks in
    %   marked4ed, marked4e and markedP4e, e.g. for
    %   MESH = MESH.refine(OBJ.marked4ed) in 2d or
    %   MESH = MESH.refine(OBJ.marked4e) for a Mesh3d.
    %
    %   MARK functions:
    %       initializeMark - Prepare empty marking arrays for a mesh
    %       mark           - Mark entities according to the strategy
    %
    %   MARK properties:
    %       initialPredict - Initial value of the predictor
    %       theta          - Bulk parameter
    %       gammaH         - Update parameter for h-refined elements
    %       gammaP         - Update parameter for p-refined elements
    %       gammaN         - Update parameter for elements not refined
    %       marked4ed      - Marks of edges to be refined
    %       marked4e       - Marks of elements to be refined
    %       markedP4e      - Marks of elements to be p-refined
    %       predict4e      - Predictor of squared indicators per element
    %
    %   See also MarkBulk, MarkMaximum, MarkMean, MarkFixedFraction,
    %   MarkHP, MarkP, MarkUniform

    % Copyright (c) 2026 by the fem4m authors

    properties
        initialPredict = inf; % initial value of the predictor
        theta = 0.5; % bulk parameter
        gammaH = 4; % update parameter for h-refined elements
        gammaP = 0.4; % update parameter for p-refined elements
        gammaN = 1; % update parameter for elements not refined
        marked4ed  % marks of edges to be refined
        marked4e   % marks of elements to be refined
        markedP4e  % marks of elements to be p-refined
        predict4e = [] % predictor of squared indicators per element
    end

    methods
        function obj = Mark(theta,gammaH,gammaP,gammaN,initialPredict)
            %MARK - Create a marking strategy object
            if ~(nargin < 1 || isempty(theta))
                obj.theta = theta;
            end
            if ~(nargin < 2 || isempty(gammaH))
                obj.gammaH = gammaH;
            end
            if ~(nargin < 3 || isempty(gammaP))
                obj.gammaP = gammaP;
            end
            if ~(nargin < 4 || isempty(gammaN))
                obj.gammaN = gammaN;
            end
            if ~(nargin < 5 || isempty(initialPredict))
                obj.initialPredict = initialPredict;
            end
        end

        function obj = initializeMark(obj,mesh)
            %initializeMark - Prepare empty marking arrays for a mesh
            %   OBJ = initializeMark(OBJ,MESH) resets marked4ed, marked4e
            %   and markedP4e to zeros on MESH. An empty predict4e is
            %   filled with OBJ.initialPredict. When the number of
            %   elements has changed, every element takes the predictor
            %   of its parent MESH.parent4e, or OBJ.initialPredict when
            %   MESH has no matching parent4e (e.g. MeshPoly classes or
            %   after coarsening).
            obj.marked4ed = zeros(mesh.nrEdges,1);
            obj.marked4e = zeros(mesh.nrElems,1);
            obj.markedP4e = zeros(mesh.nrElems,1);
            if isempty(obj.predict4e)
                obj.predict4e = obj.initialPredict.*ones(mesh.nrElems,1);
            elseif numel(obj.predict4e) ~= mesh.nrElems
                if isprop(mesh,'parent4e') && numel(mesh.parent4e) == mesh.nrElems ...
                        && max(mesh.parent4e) <= numel(obj.predict4e)
                    obj.predict4e = obj.predict4e(mesh.parent4e);
                else
                    % no parents (e.g. after coarsening): start again
                    obj.predict4e = obj.initialPredict.*ones(mesh.nrElems,1);
                end
            end
        end
    end

    methods (Abstract)
        %MARK - Mark entities according to the strategy
        obj = mark(obj)
    end
end
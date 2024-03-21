classdef MarkHP < Mark
    %MarkHP - NVB based hp marking strategy bisecting one edge
    %   OBJ = MarkHP() creates an hp marking strategy with the default
    %   parameters of Mark.
    %
    %   OBJ = MarkHP(THETA,GAMMAH) sets the update parameter GAMMAH
    %   of the predictor for h-refined elements as described in Mark;
    %   THETA is ignored.
    %
    %   OBJ = MarkHP(...,GAMMAP) also sets the update parameter GAMMAP
    %   for p-refined elements.
    %
    %   OBJ = MarkHP(...,GAMMAN) also sets the update parameter GAMMAN
    %   for elements that are not refined.
    %
    %   OBJ = MarkHP(...,INITIALPREDICT) also sets the initial value
    %   INITIALPREDICT of the predictor.
    %
    %   MarkHP functions:
    %       mark - Mark elements for h- or p-refinement by NVB
    %
    %   See also Mark, MarkHP3Edges, MarkP

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = mark(obj,fe,varargin)
            %MARK - Mark elements for h- or p-refinement by NVB
            %   OBJ = MARK(OBJ,FE,ETA4E) selects the elements whose
            %   squared indicator ETA4E.^2 is at least three quarters of
            %   the mean squared indicator. A selected element with
            %   ETA4E.^2 greater than its predictor OBJ.predict4e is
            %   marked for h-refinement in OBJ.marked4e, together with its
            %   refinement edge (the first edge of FE.mesh.ed4e) in
            %   OBJ.marked4ed for newest vertex bisection; the other
            %   selected elements are marked for p-refinement in
            %   OBJ.markedP4e. The predictor becomes
            %   OBJ.gammaH/2*(1/2)^P*ETA4E.^2 on the h-refined elements,
            %   P = FE.degree4e the elementwise polynomial degree (e.g.
            %   of HPk), OBJ.gammaP*ETA4E.^2 on the p-refined elements and
            %   is multiplied by OBJ.gammaN on the elements not selected.
            %   ETA4E is a column vector of indicators per element of a
            %   triangular mesh. Refine the mesh with
            %   MESH = MESH.refine('nvb',OBJ.marked4ed).
            %
            %   OBJ = MARK(OBJ,FE,ETA4E,INDICATORS2,...) ignores all
            %   indicator vectors after ETA4E.
            mesh = fe.mesh;
            obj = obj.initializeMark(mesh);
            indicators = varargin{1};

            I = indicators.^2 >= 0.75*sum(indicators.^2)/mesh.nrElems;

            J =  I & indicators.^2 > obj.predict4e;
            K =  I & ~(indicators.^2 > obj.predict4e);

            obj.predict4e(~I) = obj.gammaN*obj.predict4e(~I);
            obj.predict4e(J) = 1/2*obj.gammaH*(1/2).^(fe.degree4e(J)).*indicators(J).^2;
            obj.predict4e(K) = obj.gammaP*indicators(K).^2;

            obj.marked4ed(mesh.ed4e(J,1)) = 1;
            obj.marked4e(J) = 1;
            obj.markedP4e(K) = 1;
        end
    end
end
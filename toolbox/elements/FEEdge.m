classdef FEEdge
    %FEEdge - Abstract base class of edge-based finite element objects
    %   FEEdge defines the interface of finite element objects whose
    %   degrees of freedom live on the edges of the mesh (dof4ed) and
    %   evaluates their finite element functions from the edge basis
    %   getEdgeBasis(OBJ,PTS,X,PARTS) of the subclass.
    %
    %   FEEdge functions:
    %       getEdgeValue     - Evaluate the edge finite element function
    %       getEdgeGradValue - Evaluate the gradient of the edge function
    %
    %   FEEdge properties:
    %       dof4ed - Local-to-global degrees of freedom per edge
    %
    %   See also FE, FEEdgeVector, FEFace

    % Copyright (c) 2026 by the fem4m authors

    properties(Abstract)
        dof4ed    % local-to-global degrees of freedom per edge
    end

    methods

        function val = getEdgeValue(obj,~,x,parts,c)
            %getEdgeValue - Evaluate the edge finite element function
            %   VAL = getEdgeValue(OBJ,PTS,X,PARTS,C) evaluates the edge
            %   finite element function with coefficients C on the edges
            %   PARTS at the point X (1-by-1) of the reference edge [0,1].
            %   PTS is ignored. VAL has one row per edge and the values on
            %   the two sides of the edges along the last dimension
            %   (numel(PARTS)-by-2 for a scalar edge element).
            val = obj.getEdgeBasis([],x,parts); %#ok<MCNPN>
            c4dof = c(obj.dof4ed(parts,:));
            val = permute(sum(val.*c4dof,2),[1 3 4 2]);
        end

        function val = getEdgeGradValue(obj,~,x,parts,c)
            %getEdgeGradValue - Evaluate the gradient of the edge function
            %   VAL = getEdgeGradValue(OBJ,PTS,X,PARTS,C) evaluates the
            %   gradient of the edge finite element function with
            %   coefficients C on the edges PARTS at the point X (1-by-1)
            %   of the reference edge [0,1] from the edge gradients
            %   getEdgeGradBasis, which the subclass must provide. PTS is
            %   ignored.
            val = obj.getEdgeGradBasis([],x,parts); %#ok<MCNPN>
            c4dof = c(obj.dof4ed(parts,:));
            val = permute(sum(val.*c4dof,2),[1 3 4 2]);
        end

    end
end
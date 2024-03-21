classdef FEFace
    %FEFace - Abstract base class of face-based finite element objects
    %   FEFace defines the interface of finite element objects whose
    %   degrees of freedom live on the faces of the mesh (dof4f) and
    %   evaluates their finite element functions from the face basis
    %   getFaceBasis(OBJ,PTS,PTS_REF,PARTS) of the subclass.
    %
    %   FEFace functions:
    %       getFaceValue - Evaluate the face finite element function
    %
    %   FEFace properties:
    %       dof4f - Local-to-global degrees of freedom per face
    %
    %   See also FE, FEEdge, DGDkFace

    % Copyright (c) 2026 by the fem4m authors

    properties(Abstract)
        dof4f    % local-to-global degrees of freedom per face
    end

    methods

        function val = getFaceValue(obj,pts,pts_ref,parts,c)
            %getFaceValue - Evaluate the face finite element function
            %   VAL = getFaceValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the
            %   face finite element function with coefficients C on the
            %   faces PARTS at the point PTS_REF (1-by-2) of the reference
            %   triangle; PTS is passed on to getFaceBasis. VAL has one
            %   row per face and one column per side of the face.
            dof4e = obj.dof4f(parts,:);
            coeff4e = c(dof4e);
            basis4e = obj.getFaceBasis(pts,pts_ref,parts); %#ok<MCNPN>
            val = permute(sum(basis4e.*coeff4e,2),[1 3 2]);
        end

    end
end
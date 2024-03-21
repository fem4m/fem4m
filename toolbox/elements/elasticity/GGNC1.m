classdef GGNC1 < AWNC1
    %GGNC1 - Gopalakrishnan-Guzman non-conforming stress element
    %   OBJ = GGNC1(DEGREE,MESH) creates the quadratic
    %   Gopalakrishnan-Guzman non-conforming stress finite element for
    %   elasticity on the triangular mesh MESH, whose space consists of
    %   all symmetric quadratic matrix fields. The 18 degrees of freedom
    %   per element are the 12 edge moments of AWNC1 and 6 element
    %   degrees of freedom: the mean values of the three stress
    %   components and, at each vertex, the value n_i*sigma*n_j for the
    %   unit normals n_i, n_j of its two edges. DEGREE must be 2.
    %
    %   GGNC1 functions:
    %       prepareCoefficients - Compute the basis expansion coefficients
    %
    %   See also AWNC1, AWNC0, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = GGNC1(degree,mesh)
            %GGNC1 - Create a Gopalakrishnan-Guzman element object
            if nargin < 2
                [degree,mesh] = deal([]);
            end
            obj = obj@AWNC1(degree,mesh);
        end

        %% Calculation of the Gopalakrishnan & Guzman NCFE transformation
        function obj = prepareCoefficients(obj)
            %prepareCoefficients - Compute the basis expansion coefficients
            %   OBJ = prepareCoefficients(OBJ) solves, on every element,
            %   the 18x18 Vandermonde system given by the 12 edge
            %   moments, the 3 element moments and the 3 vertex values
            %   n_i*sigma*n_j (all 18 rows are degrees of freedom, there
            %   are no constraints), and stores the expansion of the
            %   Gopalakrishnan-Guzman basis in the P2 basis of basisP2 in
            %   basisCoefficients.
            normals4e = obj.mesh.normal4e;
            nrElems = obj.mesh.nrElems;
            ed4e = obj.mesh.ed4e;
            normal4ed = obj.mesh.normal4ed;
            % Initialisation
            obj.basisCoefficients = zeros(nrElems,18,18);
            C = zeros(18,18);
            % Element degrees of freedom
            I = [eye(3,3),eye(3,3),eye(3,3)];
            C(13:15,:) = [20*I,5*I];
            % Element specific transformation
            for curElem = 1 : nrElems
                % Edge degrees of freedom
                curNormal = normal4ed(ed4e(curElem,:),:);
                N1 = [diag(curNormal(1,:)),curNormal(1,[2 1])'];
                N2 = [diag(curNormal(2,:)),curNormal(2,[2 1])'];
                N3 = [diag(curNormal(3,:)),curNormal(3,[2 1])'];
                curNormal = permute(normals4e(curElem,:,:),[2 3 1]);
                M1 = [diag(curNormal(1,:)),curNormal(1,[2 1])'];
                M2 = [diag(curNormal(2,:)),curNormal(2,[2 1])'];
                M3 = [diag(curNormal(3,:)),curNormal(3,[2 1])'];
                O = zeros(2,3);
                R = [30*N1,30*N1,O;-5*M1,5*M1,O;
                    O,30*N2,30*N2;O,-5*M2,5*M2;
                    30*N3,O,30*N3;5*M3,O,-5*M3];
                S = [10*N1,O,O;O,O,O;
                    O,10*N2,O;O,O,O;
                    O,O,10*N3;O,O,O];
                C(1:12,:) = [R,S];
                % vertex dofs n1*sigma*n3, n2*sigma*n1, n3*sigma*n2
                N1 = [curNormal(1,:).*curNormal(3,:),sum(curNormal(1,:).*curNormal(3,[2 1]),2)];
                N2 = [curNormal(2,:).*curNormal(1,:),sum(curNormal(2,:).*curNormal(1,[2 1]),2)];
                N3 = [curNormal(3,:).*curNormal(2,:),sum(curNormal(3,:).*curNormal(2,[2 1]),2)];
                O = zeros(size(N1));
                C(16:18,1:9) = 60*[N1, O, O; O, N2, O; O, O, N3];
                % Assembly
                x = (C/60)\eye(18);
                obj.basisCoefficients(curElem,:,:) = x';
            end
        end


    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of no nodal, 4 edge moments, and
            %   6 element dofs (3 element moments and 3 vertex values).
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 4;
            obj.nrBasis4e = 6;
        end

    end
end
classdef AW0 < AW1
    %AW0 - Reduced Arnold-Winther element for elasticity
    %   OBJ = AW0(DEGREE,MESH) creates the reduced cubic Arnold-Winther
    %   symmetric H(div)-conforming stress finite element on the
    %   triangular mesh MESH with 21 degrees of freedom per element:
    %   the values of sigma at the vertices (3 per vertex) and the
    %   moments of degree 0 and 1 of sigma*n on every edge (4 per
    %   edge), no element moments. The space consists of the
    %   symmetric cubic matrix fields whose divergence is a rigid
    %   motion. DEGREE must be 3.
    %
    %   AW0 functions:
    %       prepareCoefficients - Compute the basis expansion coefficients
    %       hessianP3           - Evaluate cubic basis Hessians at the centroid
    %
    %   See also AW1, AWk, HZ1, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = AW0(degree,mesh)
            %AW0 - Create a reduced Arnold-Winther element object
            if nargin < 2
                [degree,mesh] = deal([]);
            end
            obj = obj@AW1(degree,mesh);
        end

        %% Calculation of the Arnold Winther transformation
        function obj = prepareCoefficients(obj)
            %prepareCoefficients - Compute the basis expansion coefficients
            %   OBJ = prepareCoefficients(OBJ) solves the local 30x30
            %   Vandermonde system given by nodal values, edge moments,
            %   and the divergence constraints div(sigma) in P1 and
            %   eps(div(sigma)) = 0 (div(sigma) a rigid motion) to obtain
            %   the expansion of the reduced Arnold-Winther basis in the
            %   P3 basis.
            normals4e = obj.mesh.normal4e;
            nrElems = obj.mesh.nrElems;
            gradP14e = obj.gradP1;
            gradP14e = permute(matMul(obj.mesh.transformation,permute(gradP14e,[3 2 1])),[1 3 2]);
            ed4e = obj.mesh.ed4e;
            normal4ed = obj.mesh.normal4ed;
            length4ed = obj.mesh.length4ed;
            % Initialisation
            obj.basisCoefficients = zeros(nrElems,21,30);
            C = zeros(30,30);
            % Nodal degrees of freedom
            C(1:9,1:9) = 60*eye(9,9);
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
                S1 = [10*N1,O,O;O,O,O;
                    O,10*N2,O;O,O,O;
                    O,O,10*N3;O,O,O];
                S2 = [O,O,O;-M1,O,O;
                    O,O,O;O,-M2,O;
                    O,O,O;O,O,-M3];
                C(10:21,1:27) = [R,S1,S2];
                % Constraint div(sigma) in P1
                gradP1X = gradP14e(curElem,:,1);
                gradP1Y = gradP14e(curElem,:,2);
                dxW = [3*(gradP1X(1)-gradP1X(2)), -gradP1X(3), gradP1X(3), gradP1X(3);
                    gradP1X(1), 3*(gradP1X(2)-gradP1X(3)), -gradP1X(1), gradP1X(1);
                    -gradP1X(2), gradP1X(2), 3*(gradP1X(3)-gradP1X(1)), gradP1X(2)];
                dyW = [3*(gradP1Y(1)-gradP1Y(2)), -gradP1Y(3), gradP1Y(3), gradP1Y(3);
                    gradP1Y(1), 3*(gradP1Y(2)-gradP1Y(3)), -gradP1Y(1), gradP1Y(1);
                    -gradP1Y(2), gradP1Y(2), 3*(gradP1Y(3)-gradP1Y(1)), gradP1Y(2)];
                O = zeros(3,1);
                Q = [dxW(:,1),O,dyW(:,1),dxW(:,2),O,dyW(:,2),dxW(:,3),O,dyW(:,3),dxW(:,4),O,dyW(:,4);
                    O,dyW(:,1),dxW(:,1),O,dyW(:,2),dxW(:,2),O,dyW(:,3),dxW(:,3),O,dyW(:,4),dxW(:,4)];
                C(25:30,19:30) = 60*max(length4ed(ed4e(curElem,:)))*Q;
                % Constraints for reduced AW: the linear div(sigma) is a
                % rigid motion, eps(div(sigma)) = 0, which is constant and
                % taken at the centroid, [a_xx + c_xy; b_yy + c_xy;
                % a_xy + b_xy + c_xx + c_yy] for sigma = [a c; c b]
                H = obj.hessianP3(reshape(gradP14e(curElem,:,:),3,2));
                O = zeros(10,1);
                Q2 = [reshape([H(:,1),O,H(:,3)]',1,[]);
                    reshape([O,H(:,2),H(:,3)]',1,[]);
                    reshape([H(:,3),H(:,3),H(:,1)+H(:,2)]',1,[])];
                C(22:24,:) = 60*max(length4ed(ed4e(curElem,:)))^2*Q2;
                % Assembly
                x = (C/60)\[eye(21);zeros(9,21)];
                obj.basisCoefficients(curElem,:,:) = x';
            end
        end

        function val = hessianP3(~,grad)
            %hessianP3 - Evaluate cubic basis Hessians at the centroid
            %   VAL = hessianP3(OBJ,GRAD) returns the Hessians [d_xx,
            %   d_yy, d_xy] (columns) of the ten cubic basis functions of
            %   basisP3 (rows) at the centroid of a triangle whose
            %   barycentric coordinates have the gradients GRAD (3-by-2,
            %   one row per coordinate). VAL is 10-by-3; the rows of the
            %   three linear functions are zero.
            % d^2(lambda_i*lambda_j) = S(i,j) = g_i*g_j' + g_j*g_i', and a
            % cubic product gives the mean of the S of its three pairs
            S = @(i,j)[2*grad(i,1)*grad(j,1), 2*grad(i,2)*grad(j,2), ...
                grad(i,1)*grad(j,2) + grad(j,1)*grad(i,2)];
            val = zeros(10,3);
            val(4,:) = S(1,2);
            val(5,:) = S(2,3);
            val(6,:) = S(3,1);
            val(7,:) = (S(1,1) - S(2,2))/3;
            val(8,:) = (S(2,2) - S(3,3))/3;
            val(9,:) = (S(3,3) - S(1,1))/3;
            val(10,:) = (S(1,2) + S(2,3) + S(3,1))/3;
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of the 3 nodal stress values and
            %   4 edge moments.
            obj.nrBasis4n = 3;
            obj.nrBasis4ed = 4;
            obj.nrBasis4e = 0;
        end

    end
end
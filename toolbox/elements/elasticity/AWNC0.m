classdef AWNC0 < AWNC1
    %AWNC0 - Reduced non-conforming Arnold-Winther element for elasticity
    %   OBJ = AWNC0(DEGREE,MESH) creates the reduced non-conforming
    %   quadratic Arnold-Winther-type stress finite element on the
    %   triangular mesh MESH: the quadratic symmetric fields of AWNC1
    %   whose divergence is a rigid motion. Its 12 degrees of freedom
    %   per element are the zeroth and first order edge moments of
    %   AWNC1; there are no nodal and no element degrees of freedom.
    %   DEGREE must be 2.
    %
    %   AWNC0 functions:
    %       prepareCoefficients - Compute the basis expansion coefficients
    %       hessianP2           - Evaluate the Hessians of the quadratic basis
    %
    %   See also AWNC1, GGNC1, AW0, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = AWNC0(degree,mesh)
            %AWNC0 - Create a reduced non-conforming Arnold-Winther element object
            if nargin < 2
                [degree,mesh] = deal([]);
            end
            obj = obj@AWNC1(degree,mesh);
        end

        %% Calculation of the Arnold Winther transformation
        function obj = prepareCoefficients(obj)
            %prepareCoefficients - Compute the basis expansion coefficients
            %   OBJ = prepareCoefficients(OBJ) solves, on every element,
            %   the 18x18 Vandermonde system given by the 12 edge
            %   moments, the three normal-normal constraints of AWNC1 and
            %   the three constraints eps(div(sigma)) = 0 (div(sigma) a
            %   rigid motion), and stores the expansion of the reduced
            %   basis in the P2 basis of basisP2 in basisCoefficients.
            normals4e = obj.mesh.normal4e;
            nrElems = obj.mesh.nrElems;
            gradP14e = obj.gradP1;
            gradP14e = permute(matMul(obj.mesh.transformation,permute(gradP14e,[3 2 1])),[1 3 2]);
            ed4e = obj.mesh.ed4e;
            normal4ed = obj.mesh.normal4ed;
            length4ed = obj.mesh.length4ed;
            % Initialisation
            obj.basisCoefficients = zeros(nrElems,12,18);
            C = zeros(18,18);
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
                % Constraints n*sigma*n in P1 on E1,E2,E3
                N1 = [curNormal(1,:).^2,2*curNormal(1,1)*curNormal(1,2)];
                N2 = [curNormal(2,:).^2,2*curNormal(2,1)*curNormal(2,2)];
                N3 = [curNormal(3,:).^2,2*curNormal(3,1)*curNormal(3,2)];
                O = zeros(size(N1));
                C(13:15,10:18) = 60*[N1, O, O; O, N2, O; O, O, N3];
                % Constraints for reduced AW: the linear div(sigma) is a
                % rigid motion, eps(div(sigma)) = 0, which is constant,
                % [a_xx + c_xy; b_yy + c_xy; a_xy + b_xy + c_xx + c_yy]
                % for sigma = [a c; c b]
                H = obj.hessianP2(reshape(gradP14e(curElem,:,:),3,2));
                O = zeros(6,1);
                Q2 = [reshape([H(:,1),O,H(:,3)]',1,[]);
                    reshape([O,H(:,2),H(:,3)]',1,[]);
                    reshape([H(:,3),H(:,3),H(:,1)+H(:,2)]',1,[])];
                C(16:18,:) = 60*max(length4ed(ed4e(curElem,:)))^2*Q2;
                % Assembly
                x = (C/60)\[eye(12);zeros(6,12)];
                obj.basisCoefficients(curElem,:,:) = x';
            end
        end

        function val = hessianP2(~,grad)
            %hessianP2 - Evaluate the Hessians of the quadratic basis
            %   VAL = hessianP2(OBJ,GRAD) returns the constant Hessians
            %   [d_xx, d_yy, d_xy] (columns) of the six quadratic basis
            %   functions of basisP2 (rows) as a 6-by-3 array on a
            %   triangle whose barycentric coordinates have the gradients
            %   GRAD (3-by-2, one row each); the rows of the three linear
            %   functions are zero. OBJ is ignored.
            % d^2(lambda_i*lambda_j) = g_i*g_j' + g_j*g_i'
            S = @(i,j)[2*grad(i,1)*grad(j,1), 2*grad(i,2)*grad(j,2), ...
                grad(i,1)*grad(j,2) + grad(j,1)*grad(i,2)];
            val = zeros(6,3);
            val(4,:) = S(1,2);
            val(5,:) = S(2,3);
            val(6,:) = S(3,1);
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of no nodal, 4 edge moments, and
            %   no element moments.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 4;
            obj.nrBasis4e = 0;
        end

    end
end
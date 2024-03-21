classdef FormValue
    %FormValue - Finite element function operators for form assembly
    %   FormValue provides static methods that evaluate differential
    %   operators of finite element functions with given coefficients
    %   on the elements (suffix _T) and on the sides (suffix _S), i.e.
    %   the edges in 2d and the faces in 3d (see Form.sideOf),
    %   including jump and average variants. EstimatorForm and
    %   ErrorForm select them by name from the parsed form string;
    %   BilinearForm and LinearForm use the basis operators of FormBasis
    %   (LinearForm also uses typeOf).
    %
    %   The operators take the arguments of an integrand of integrate:
    %   the points PTS (one row per element or side of PARTS), the
    %   point PTS_REF of the reference element or side (see integrate)
    %   and the element or side indices PARTS, followed by the
    %   coefficient vector C (X for the elasticity operators) of the
    %   finite element function and the finite element space FE; VAL has
    %   one row per element or side of PARTS. The side operators evaluate
    %   the function on the neighboring elements of the sides; the
    %   missing second element of a boundary side has the values zero, so
    %   that jumps and means are the values of the first element there. A
    %   coefficient ALPHA or the elasticity tensor C may be empty for the
    %   identity; the coefficient MU of jumpnxcurlwxn_S is required.
    %
    %   FormValue functions:
    %       w_T              - Function values on the element
    %       wn_T             - Normal component of the function
    %       gradw_T          - Gradient of the function
    %       skeww_T          - Skew part of the function
    %       divw_T           - Divergence of the function
    %       divgradw_T       - Divergence of ALPHA times the gradient
    %       laplacew_T       - Laplacian of the function
    %       divepsw_T        - Divergence of C applied to eps(u)
    %       biharmonicw_T    - Biharmonic of the function
    %       epsw_T           - Symmetric gradient of the function
    %       symw_T           - Symmetric part of the function
    %       curlw_T          - Curl of the function
    %       rotw_T           - Rot of the function
    %       d2w_T            - Hessian of the function
    %       curlcurlw_T      - Curl-curl of the function
    %       curlrotw_T       - Curl of the rot of the function
    %       rotrotw_T        - Rot of the rot of the function
    %       symcurlw_T       - Symmetric part of the curl of the function
    %       divdivw_T        - Double divergence of the function
    %       divd2w_T         - Divergence of the Hessian of the function
    %       trw_T            - Trace of the function
    %       w_S              - Function values on the sides
    %       gradw_S          - Gradient of the function on the sides
    %       wn_S             - Normal component of the function on the sides
    %       wt_S             - Tangential component of the function
    %       wtt_S            - Tangential-tangential value of the function
    %       gradwn_S         - Normal derivative of the function
    %       d2wn_S           - Hessian of the function times the normal
    %       divwn_S          - Divergence of the function times the normal
    %       divd2wn_S        - Normal derivative of the Laplacian
    %       epswn_S          - Normal symmetric gradient
    %       rotwt_S          - Tangential rot of the function
    %       gradwntt_S       - Tangential derivative of n'*w*t
    %       jumpw_S          - Jump of the function across the sides
    %       jumpgradw_S      - Jump of the function gradient
    %       jumpwn_S         - Normal jump of the function
    %       jumpwt_S         - Tangential jump of the function
    %       jumpgradwn_S     - Normal jump of the function gradient
    %       jumpwnn_S        - Second normal jump of the function
    %       jumpwtn_S        - Tangential-normal jump of the function
    %       jumpwnt_S        - Normal-tangential jump of the function
    %       jumpwtt_S        - Tangential-tangential jump of the function
    %       jumpepswn_S      - Jump of the symmetric gradient
    %       jumpd2wnn_S      - Second normal jump of the Hessian
    %       jumpd2wtn_S      - Mixed second jump of the Hessian
    %       jumprotwt_S      - Tangential jump of the rot
    %       jumpdivd2wn_S    - Jump of the normal derivative of the Laplacian
    %       jumpdivd2wt_S    - Jump of the tangential derivative of laplace(u)
    %       jumpd2wt_S       - Tangential jump of the Hessian
    %       jumpsymcurlwt_S  - Tangential jump of the symmetric curl
    %       jumpdivw_S       - Jump of the divergence
    %       jumpdivwn_S      - Normal jump of the divergence
    %       jumpdivwt_S      - Tangential jump of the divergence
    %       jumpd3wnnn_S     - Third normal jump of the function
    %       jumpgradwntt_S   - Jump of the tangential derivative of n'*w*t
    %       jumpwxn_S        - Jump times the side normals
    %       jumpnxcurlwxn_S  - Jump of n x (MU*curl(u)) x n
    %       averagecurlw_S   - Mean of the curl across the sides
    %       averagewtt_S     - Mean of the tangential-tangential value
    %       averagegradwtt_S - Mean of the tangential gradient
    %       averaged2wt_S    - Mean of the Hessian times the tangent
    %       averagedivd2wn_S - Mean normal derivative of the Laplacian
    %       averagedivd2wt_S - Mean tangential derivative of the Laplacian
    %       averagedivwn_S   - Mean of the divergence times the normal
    %       averagedivwt_S   - Mean of the divergence times the tangent
    %       sideValues       - Function values on both neighboring elements
    %       oneSide          - Side operator on the first neighboring element
    %       jump             - Jump of a side operator across the sides
    %       average          - Mean of a side operator across the sides
    %       typeOf           - Value type of function values
    %       oneSideTimes     - One-sided values times directions
    %       contractOneSide  - One-sided matrix values between two directions
    %       getJump          - Difference of the two neighboring values
    %       getAverage       - Mean of the two neighboring values
    %
    %   See also FormBasis, Form, EstimatorForm, ErrorForm, FormEval,
    %   FormType

    % Copyright (c) 2026 by the fem4m authors

    methods(Static)

        %% on volume
        function val = w_T(pts,pts_ref,parts,c,fe)
            %w_T - Function values on the element
            %   VAL = w_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   finite element function with coefficients C on FE.
            val = fe.getValue(pts,pts_ref,parts,c);
        end

        function val = wn_T(pts,pts_ref,parts,c,fe)
            %wn_T - Normal component of the function
            %   VAL = wn_T(PTS,PTS_REF,PARTS,C,FE) returns the normal component
            %   of the finite element function by FE.getNormalValue, for elements
            %   that provide getNormalValue (no element of the toolbox does).
            val = fe.getNormalValue(pts,pts_ref,parts,c);
        end

        function val = gradw_T(pts,pts_ref,parts,c,fe)
            %gradw_T - Gradient of the function
            %   VAL = gradw_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   gradient of the finite element function.
            val = fe.getGradValue(pts,pts_ref,parts,c);
        end

        function val = skeww_T(pts,pts_ref,parts,c,fe)
            %skeww_T - Skew part of the function
            %   VAL = skeww_T(PTS,PTS_REF,PARTS,C,FE) evaluates the skew part
            %   (M-M')/2 of the matrix-valued finite element function M, stored
            %   column by column like M.
            val = fe.getValue(pts,pts_ref,parts,c);
            val = reshape(val,size(val,1),sqrt(size(val,2)),sqrt(size(val,2)));
            val = 1/2*(val-permute(val,[1 3 2]));
            val = val(:,:);
        end

        function val = divw_T(pts,pts_ref,parts,c,fe)
            %divw_T - Divergence of the function
            %   VAL = divw_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   divergence of the finite element function.
            val = fe.getDivValue(pts,pts_ref,parts,c);
        end

        function val = divgradw_T(pts,pts_ref,parts,c,fe,alpha)
            %divgradw_T - Divergence of ALPHA times the gradient
            %   VAL = divgradw_T(PTS,PTS_REF,PARTS,C,FE,ALPHA) evaluates
            %   div(ALPHA*grad(u)) of the scalar finite element function u with
            %   coefficients C from its gradient and second derivatives. ALPHA
            %   is empty (identity), a scalar or a function handle (see
            %   FormEval.evalCoefficient) for the d-by-d diffusion matrix, or a
            %   cell {KAPPA,DKAPPA} with such a coefficient or a constant d-by-d
            %   matrix KAPPA and the function handle DKAPPA of the points that
            %   returns the derivatives of KAPPA as an N-by-d-by-d-by-d array
            %   with DKAPPA(:,i,j,k) = d KAPPA_ij/dx_k. The product-rule term
            %   with the derivatives of the coefficient is only included for
            %   {KAPPA,DKAPPA}; otherwise VAL is the sum of ALPHA_ij*d_i d_j u.
            gradVal = fe.getGradValue(pts,pts_ref,parts,c);
            d2Val = fe.getD2Value(pts,pts_ref,parts,c);
            if isempty(alpha) || isscalar(alpha)
                alpha4e = FormEval.evalCoefficient(pts,pts_ref,parts,alpha,FormType.vector(size(pts,2)));
            else
                alpha4e = FormEval.evalCoefficient(pts,pts_ref,parts,alpha{1},FormType.vector(size(pts,2)));
            end
            if ~isempty(alpha) && numel(alpha)==2
                dalpha4e = alpha{2};
                dalpha4e = dalpha4e(pts);
            else
                dalpha4e = zeros(size(pts,1),size(pts,2),size(pts,2),size(pts,2));
            end
            % div(kappa gradU)
            if size(pts,2) == 2
                val = ...
                    alpha4e(:,1,1).*d2Val(:,1) + alpha4e(:,1,2).*d2Val(:,3) ...
                    + alpha4e(:,2,1).*d2Val(:,3) + alpha4e(:,2,2).*d2Val(:,2) ...
                    + dalpha4e(:,1,1,1).*gradVal(:,1) + dalpha4e(:,1,2,1).*gradVal(:,2) ...
                    + dalpha4e(:,2,1,2).*gradVal(:,1) + dalpha4e(:,2,2,2).*gradVal(:,2);
            else
                val = ...
                    alpha4e(:,1,1).*d2Val(:,1) + alpha4e(:,1,2).*d2Val(:,4) + alpha4e(:,1,3).*d2Val(:,5) ...
                    +alpha4e(:,2,1).*d2Val(:,4) + alpha4e(:,2,2).*d2Val(:,2) + alpha4e(:,2,3).*d2Val(:,6) ...
                    +alpha4e(:,3,1).*d2Val(:,5) + alpha4e(:,3,2).*d2Val(:,6) + alpha4e(:,3,3).*d2Val(:,3) ...
                    + dalpha4e(:,1,1,1).*gradVal(:,1) + dalpha4e(:,1,2,1).*gradVal(:,2) +  dalpha4e(:,1,3,1).*gradVal(:,3) ...
                    + dalpha4e(:,2,1,2).*gradVal(:,1) + dalpha4e(:,2,2,2).*gradVal(:,2) +  dalpha4e(:,2,3,2).*gradVal(:,3) ...
                    + dalpha4e(:,3,1,3).*gradVal(:,1) + dalpha4e(:,3,2,3).*gradVal(:,2) +  dalpha4e(:,3,3,3).*gradVal(:,3);
            end
        end

        function val = laplacew_T(pts,pts_ref,parts,c,fe)
            %laplacew_T - Laplacian of the function
            %   VAL = laplacew_T(PTS,PTS_REF,PARTS,C,FE) evaluates
            %   the Laplacian of the finite element function.
            val = fe.getLaplaceValue(pts,pts_ref,parts,c);
        end

        function val = divepsw_T(pts,pts_ref,parts,x,fe,C)
            %divepsw_T - Divergence of C applied to eps(u)
            %   VAL = divepsw_T(PTS,PTS_REF,PARTS,X,FE,C) evaluates
            %   div(C*eps(u)) of the vector-valued finite element function u
            %   with coefficients X from its second derivatives. C is a
            %   coefficient (see FormEval.evalCoefficient) for the 3-by-3 (2d)
            %   or 6-by-6 (3d) matrix that acts on the symmetric gradient stored
            %   as [11 22 12] or [11 22 33 12 13 23]. Only lambda = C(1,2) and
            %   mu = C(end,end)/2 are read, i.e. C is taken as isotropic and VAL
            %   is mu*laplace(u) + (lambda+mu)*grad(div(u)); an empty C (the
            %   identity) gives div(eps(u)).
            if isa(fe.mesh,'Mesh3d')
                C = FormEval.evalCoefficient(pts,pts_ref,parts,C,FormType.symmatrix3d);
                lambda = C(:,1,2);
                mu = C(:,6,6)./2;
                d2 = fe.getD2Value(pts,pts_ref,parts,x);
                val = [lambda*(d2(:,2,4)+d2(:,3,5))+(lambda+2*mu)*d2(:,1,1)+mu*(d2(:,1,2)+d2(:,2,4)+d2(:,1,3)+d2(:,3,5)),...
                    lambda*(d2(:,1,4)+d2(:,3,6))+(lambda+2*mu)*d2(:,2,2)+mu*(d2(:,1,4)+d2(:,2,1)+d2(:,2,3)+d2(:,3,6)),...
                    lambda*(d2(:,1,5)+d2(:,2,6))+(lambda+2*mu)*d2(:,3,3)+mu*(d2(:,1,5)+d2(:,3,1)+d2(:,2,6)+d2(:,3,2))];
            else
                C = FormEval.evalCoefficient(pts,pts_ref,parts,C,FormType.symmatrix2d);
                lambda = C(:,1,2);
                mu = C(:,3,3)./2;
                d2 = fe.getD2Value(pts,pts_ref,parts,x);
                val = [d2(:,1,1).*(lambda+2*mu)+d2(:,1,2).*mu+d2(:,2,3).*(lambda+mu),...
                    d2(:,2,2).*(lambda+2*mu)+d2(:,2,1).*mu+d2(:,1,3).*(lambda+mu)];
            end
        end

        function val = biharmonicw_T(pts,pts_ref,parts,c,fe)
            %biharmonicw_T - Biharmonic of the function
            %   VAL = biharmonicw_T(PTS,PTS_REF,PARTS,C,FE)
            %   evaluates the biharmonic of the finite element
            %   function.
            val = fe.getBiharmonicValue(pts,pts_ref,parts,c);
        end

        function val = epsw_T(pts,pts_ref,parts,c,fe)
            %epsw_T - Symmetric gradient of the function
            %   VAL = epsw_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   symmetric gradient of the finite element function.
            val = fe.getSymGradValue(pts,pts_ref,parts,c);
        end

        function val = symw_T(pts,pts_ref,parts,c,fe)
            %symw_T - Symmetric part of the function
            %   VAL = symw_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   symmetric part (M+M')/2 of the matrix-valued function M,
            %   stored as [11 22 12] in 2d and [11 22 33 12 13 23] in 3d
            %   (see FormType.symPart).
            val = fe.getValue(pts,pts_ref,parts,c);
            val = FormValue.typeOf(val,fe).symPart(val,2);
        end

        function val = curlw_T(pts,pts_ref,parts,c,fe)
            %curlw_T - Curl of the function
            %   VAL = curlw_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   curl of the finite element function: the vector curl
            %   (d_y u, -d_x u) of a 2d scalar function and the curl of
            %   a 3d vector function.
            val = fe.getCurlValue(pts,pts_ref,parts,c);
        end

        function val = rotw_T(pts,pts_ref,parts,c,fe)
            %rotw_T - Rot of the function
            %   VAL = rotw_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   rotation d_x u_2 - d_y u_1 of the 2d vector-valued finite
            %   element function, row-wise for matrix-valued ones.
            val = fe.getRotValue(pts,pts_ref,parts,c);
        end

        function val = d2w_T(pts,pts_ref,parts,c,fe)
            %d2w_T - Hessian of the function
            %   VAL = d2w_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   second derivatives of the finite element function.
            val = fe.getD2Value(pts,pts_ref,parts,c);
        end

        function val = curlcurlw_T(pts,pts_ref,parts,c,fe)
            %curlcurlw_T - Curl-curl of the function
            %   VAL = curlcurlw_T(PTS,PTS_REF,PARTS,C,FE) evaluates
            %   the curl-curl of the finite element function: the Airy
            %   matrix curl(curl(u)) = [u_yy -u_xy; -u_xy u_xx] of a 2d
            %   scalar function and the curl-curl of a 3d vector
            %   function.
            val = fe.getCurlCurlValue(pts,pts_ref,parts,c);
        end

        function val = curlrotw_T(pts,pts_ref,parts,c,fe)
            %curlrotw_T - Curl of the rot of the function
            %   VAL = curlrotw_T(PTS,PTS_REF,PARTS,C,FE) evaluates
            %   curl(rot(u)) of the 2d vector-valued finite element
            %   function.
            val = fe.getCurlRotValue(pts,pts_ref,parts,c);
        end

        function val = rotrotw_T(pts,pts_ref,parts,c,fe)
            %rotrotw_T - Rot of the rot of the function
            %   VAL = rotrotw_T(PTS,PTS_REF,PARTS,C,FE) evaluates
            %   the scalar rot(rot(sigma)) of the 2d symmetric
            %   matrix-valued finite element function.
            val = fe.getRotRotValue(pts,pts_ref,parts,c);
        end

        function val = symcurlw_T(pts,pts_ref,parts,c,fe)
            %symcurlw_T - Symmetric part of the curl of the function
            %   VAL = symcurlw_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   symmetric part of the row-wise curl, with the rows
            %   (d_y v_i, -d_x v_i), of the 2d vector-valued finite element
            %   function v: [d_y v_1, -d_x v_2, (d_y v_2 - d_x v_1)/2].
            grad = fe.getGradValue(pts,pts_ref,parts,c);
            if size(fe.mesh.c4n,2) ~= 2 || ~ismatrix(grad) || size(grad,2) ~= 4
                error('fem4m:FormValue:symCurlNotDefined', ...
                    'FormValue: sym(curl(u)) is only defined for 2d vector fields');
            end
            % grad = [dx v1, dx v2, dy v1, dy v2], curl column by column
            curl = [grad(:,3:4),-grad(:,1:2)];
            val = FormType.matrix2d.symPart(curl,2);
        end

        function val = divdivw_T(pts,pts_ref,parts,c,fe)
            %divdivw_T - Double divergence of the function
            %   VAL = divdivw_T(PTS,PTS_REF,PARTS,C,FE) evaluates the
            %   double divergence sum_ij d_i d_j s_ij of the symmetric
            %   matrix-valued finite element function s.
            val = fe.getDivDivValue(pts,pts_ref,parts,c);
        end

        function val = divd2w_T(pts,pts_ref,parts,c,fe)
            %divd2w_T - Divergence of the Hessian of the function
            %   VAL = divd2w_T(PTS,PTS_REF,PARTS,C,FE) evaluates
            %   div(D^2 u) = grad(laplace(u)) of the scalar finite element
            %   function from the third derivatives stored as [xxx yyy xxy
            %   xyy] in 2d and [xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz] in
            %   3d.
            d3 = fe.getD3Value(pts,pts_ref,parts,c);
            if size(d3,2) == 4
                val = [sum(d3(:,[1 4]),2),sum(d3(:,[2 3]),2)];
            else
                val = [sum(d3(:,[1 6 8]),2),sum(d3(:,[4 2 9]),2),sum(d3(:,[5 7 3]),2)];
            end
        end

        function val = trw_T(pts,pts_ref,parts,c,fe)
            %trw_T - Trace of the function
            %   VAL = trw_T(PTS,PTS_REF,PARTS,C,FE) evaluates the trace of the
            %   matrix-valued finite element function: by FE.getTraceValue for
            %   elements that provide getTraceBasis (no element of the toolbox
            %   does), otherwise as the sum of the diagonal entries.
            if ismethod(fe,'getTraceBasis')
                val = fe.getTraceValue(pts,pts_ref,parts,c);
            else
                val = fe.getValue(pts,pts_ref,parts,c);
                type = FormValue.typeOf(val,fe);
                val = sum(val(:,type.diagonal),2);
            end
        end

        %% on sides (edges in 2d, faces in 3d)
        function val = w_S(pts,pts_ref,parts,c,fe)
            %w_S - Function values on the sides
            %   VAL = w_S(PTS,PTS_REF,PARTS,C,FE) evaluates the finite
            %   element function with coefficients C of FE on the first
            %   neighboring element of the sides PARTS.
            val = FormValue.oneSide('','',pts,pts_ref,parts,c,fe);
        end

        function val = gradw_S(pts,pts_ref,parts,c,fe,alpha)
            %gradw_S - Gradient of the function on the sides
            %   VAL = gradw_S(PTS,PTS_REF,PARTS,C,FE,ALPHA) applies the
            %   coefficient ALPHA to the gradient of the finite element
            %   function on the first neighboring element of the sides PARTS:
            %   a d-by-d matrix for scalar functions and a matrix acting on all
            %   gradient components otherwise (see oneSide).
            val = FormValue.oneSide('Grad','',pts,pts_ref,parts,c,fe,alpha);
        end

        function val = wn_S(pts,pts_ref,parts,c,fe,alpha)
            %wn_S - Normal component of the function on the sides
            %   VAL = wn_S(PTS,PTS_REF,PARTS,C,FE,ALPHA) multiplies ALPHA
            %   applied to the function values by the side normals (see
            %   oneSideTimes): u*n of scalar, the normal component of vector
            %   and the normal traction M*n of (symmetric) matrix functions.
            val = FormValue.oneSide('','n',pts,pts_ref,parts,c,fe,alpha);
        end

        function val = wt_S(pts,pts_ref,parts,c,fe,alpha)
            %wt_S - Tangential component of the function
            %   VAL = wt_S(PTS,PTS_REF,PARTS,C,FE,ALPHA) multiplies ALPHA
            %   applied to the function values by the edge tangents (see
            %   oneSideTimes): u*t of scalar, the tangential component of
            %   vector and M*t of matrix functions (edges only).
            val = FormValue.oneSide('','t',pts,pts_ref,parts,c,fe,alpha);
        end

        function val = wtt_S(pts,pts_ref,parts,c,fe)
            %wtt_S - Tangential-tangential value of the function
            %   VAL = wtt_S(PTS,PTS_REF,PARTS,C,FE) returns t'*M*t of the
            %   matrix-valued function M on the first neighboring element of the
            %   edges PARTS with the edge tangents t (edges only).
            val = FormValue.oneSide('','tt',pts,pts_ref,parts,c,fe);
        end

        function val = gradwn_S(pts,pts_ref,parts,c,fe,alpha)
            %gradwn_S - Normal derivative of the function
            %   VAL = gradwn_S(PTS,PTS_REF,PARTS,C,FE,ALPHA) multiplies ALPHA
            %   applied to the function gradient by the side normals (see
            %   oneSide): ALPHA*grad(u)*n of scalar functions and the gradient
            %   matrix times n of vector functions.
            val = FormValue.oneSide('Grad','n',pts,pts_ref,parts,c,fe,alpha);
        end

        function val = d2wn_S(pts,pts_ref,parts,c,fe)
            %d2wn_S - Hessian of the function times the normal
            %   VAL = d2wn_S(PTS,PTS_REF,PARTS,C,FE) multiplies the second
            %   derivatives of the scalar function on the first neighboring
            %   element of the sides PARTS by the normal of that element,
            %   D^2 u*n.
            val = FormValue.oneSide('D2','n',pts,pts_ref,parts,c,fe);
        end

        function val = divwn_S(pts,pts_ref,parts,c,fe)
            %divwn_S - Divergence of the function times the normal
            %   VAL = divwn_S(PTS,PTS_REF,PARTS,C,FE) multiplies the
            %   divergence of the function on the first neighboring element
            %   of the sides PARTS by the normal of that element: div(u)*n
            %   for vector and div(s)*n for matrix-valued functions.
            val = FormValue.oneSide('Div','n',pts,pts_ref,parts,c,fe);
        end

        function val = divd2wn_S(pts,pts_ref,parts,c,fe)
            %divd2wn_S - Normal derivative of the Laplacian
            %   VAL = divd2wn_S(PTS,PTS_REF,PARTS,C,FE) returns div(D^2 u)*n
            %   = grad(laplace(u))*n of the scalar function on the first
            %   neighboring element of the sides PARTS with the normal of
            %   that element.
            val = FormValue.oneSide('DivD2','n',pts,pts_ref,parts,c,fe);
        end

        function val = epswn_S(pts,pts_ref,parts,x,fe,C)
            %epswn_S - Normal symmetric gradient
            %   VAL = epswn_S(PTS,PTS_REF,PARTS,X,FE,C) applies the elasticity
            %   tensor C to the symmetric gradient of the vector function with
            %   coefficients X on the first neighboring element of the sides
            %   PARTS and multiplies the result by the side normals,
            %   (C*eps(u))*n (see oneSide).
            val = FormValue.oneSide('SymGrad','n',pts,pts_ref,parts,x,fe,C);
        end

        function val = rotwt_S(pts,pts_ref,parts,c,fe)
            %rotwt_S - Tangential rot of the function
            %   VAL = rotwt_S(PTS,PTS_REF,PARTS,C,FE) multiplies the rot of the
            %   function by the edge tangents (see oneSideTimes), e.g. rot(M)*t
            %   of the row-wise rot of a matrix function M (edges only).
            val = FormValue.oneSide('Rot','t',pts,pts_ref,parts,c,fe);
        end

        function val = gradwntt_S(pts,pts_ref,parts,c,fe)
            %gradwntt_S - Tangential derivative of n'*w*t
            %   VAL = gradwntt_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   tangential derivative d_t(n'*w*t) of the symmetric
            %   matrix-valued function w with the edge normals n and
            %   the tangents t = (-n_2,n_1) (edges of 2d meshes only).
            Form.checkSides(fe.mesh,'planarEdges','FormValue','gradwntt_S');
            val = FormValue.sideValues('Grad',pts,pts_ref,parts,c,fe);
            val = Form.getT1(val);
            normals = Form.getNormal4ed(pts_ref,parts,fe.mesh);
            normals = Form.getT1(normals);
            n1 = normals(:,1); n2 = normals(:,2);
            val = val(:,:).*[n1.*n2.^2,-n1.*n2.^2,n2.^3-n1.^2.*n2,-n1.^2.*n2,n1.^2.*n2,-n1.*n2.^2+n1.^3];
            val = sum(val,2);
        end

        %% jumps & averages (value)
        function val = jumpw_S(pts,pts_ref,parts,c,fe)
            %jumpw_S - Jump of the function across the sides
            %   VAL = jumpw_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   difference of the function values on the two
            %   neighboring elements of each side in PARTS.
            val = FormValue.jump('','',pts,pts_ref,parts,c,fe);
        end

        function val = jumpgradw_S(pts,pts_ref,parts,c,fe)
            %jumpgradw_S - Jump of the function gradient
            %   VAL = jumpgradw_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   difference of the function gradient on the two
            %   neighboring elements of each side in PARTS.
            val = FormValue.jump('Grad','',pts,pts_ref,parts,c,fe);
        end

        function val = jumpwn_S(pts,pts_ref,parts,c,fe,alpha)
            %jumpwn_S - Normal jump of the function
            %   VAL = jumpwn_S(PTS,PTS_REF,PARTS,C,FE,ALPHA) returns the
            %   jump of ALPHA applied to the function values times the
            %   side normals (see jump).
            val = FormValue.jump('','n',pts,pts_ref,parts,c,fe,alpha);
        end

        function val = jumpwt_S(pts,pts_ref,parts,c,fe,alpha)
            %jumpwt_S - Tangential jump of the function
            %   VAL = jumpwt_S(PTS,PTS_REF,PARTS,C,FE,ALPHA) returns the
            %   jump of ALPHA applied to the function values times the
            %   edge tangents (edges only, see jump).
            val = FormValue.jump('','t',pts,pts_ref,parts,c,fe,alpha);
        end

        function val = jumpgradwn_S(pts,pts_ref,parts,c,fe,alpha)
            %jumpgradwn_S - Normal jump of the function gradient
            %   VAL = jumpgradwn_S(PTS,PTS_REF,PARTS,C,FE,ALPHA) returns
            %   the jump of ALPHA applied to the function gradient times
            %   the side normals (see jump).
            val = FormValue.jump('Grad','n',pts,pts_ref,parts,c,fe,alpha);
        end

        function val = jumpwnn_S(pts,pts_ref,parts,c,fe)
            %jumpwnn_S - Second normal jump of the function
            %   VAL = jumpwnn_S(PTS,PTS_REF,PARTS,C,FE) returns the jump
            %   of the matrix-valued function tested twice against the
            %   side normals (see jump).
            val = FormValue.jump('','nn',pts,pts_ref,parts,c,fe);
        end

        function val = jumpwtn_S(pts,pts_ref,parts,c,fe)
            %jumpwtn_S - Tangential-normal jump of the function
            %   VAL = jumpwtn_S(PTS,PTS_REF,PARTS,C,FE) returns the jump
            %   of n'*M*t of the matrix-valued function M across the
            %   edges PARTS (edges only, see jump).
            val = FormValue.jump('','tn',pts,pts_ref,parts,c,fe);
        end

        function val = jumpwnt_S(pts,pts_ref,parts,c,fe)
            %jumpwnt_S - Normal-tangential jump of the function
            %   VAL = jumpwnt_S(PTS,PTS_REF,PARTS,C,FE) returns the jump
            %   of t'*M*n of the matrix-valued function M across the
            %   edges PARTS (edges only, see jump).
            val = FormValue.jump('','nt',pts,pts_ref,parts,c,fe);
        end

        function val = jumpwtt_S(pts,pts_ref,parts,c,fe)
            %jumpwtt_S - Tangential-tangential jump of the function
            %   VAL = jumpwtt_S(PTS,PTS_REF,PARTS,C,FE) returns the jump
            %   of the matrix-valued function tested twice against the
            %   edge tangents (edges only, see jump).
            val = FormValue.jump('','tt',pts,pts_ref,parts,c,fe);
        end

        function val = jumpepswn_S(pts,pts_ref,parts,x,fe,C)
            %jumpepswn_S - Jump of the symmetric gradient
            %   VAL = jumpepswn_S(PTS,PTS_REF,PARTS,X,FE,C) returns the jump
            %   of (C*eps(u))*n, the elasticity tensor C applied to the
            %   symmetric gradient of the vector function u with coefficients
            %   X times the side normals n (see jump, epswn_S).
            val = FormValue.jump('SymGrad','n',pts,pts_ref,parts,x,fe,C);
        end

        function val = jumpd2wnn_S(pts,pts_ref,parts,c,fe)
            %jumpd2wnn_S - Second normal jump of the Hessian
            %   VAL = jumpd2wnn_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   jump of the second derivative of the function, tested
            %   twice against the side normals (see jump).
            val = FormValue.jump('D2','nn',pts,pts_ref,parts,c,fe);
        end

        function val = jumpd2wtn_S(pts,pts_ref,parts,c,fe)
            %jumpd2wtn_S - Mixed second jump of the Hessian
            %   VAL = jumpd2wtn_S(PTS,PTS_REF,PARTS,C,FE) returns the jump of
            %   n'*D^2 u*t of the scalar function with the edge tangents t and
            %   normals n (edges only, see jump).
            val = FormValue.jump('D2','tn',pts,pts_ref,parts,c,fe);
        end

        function val = jumprotwt_S(pts,pts_ref,parts,c,fe)
            %jumprotwt_S - Tangential jump of the rot
            %   VAL = jumprotwt_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   jump of the rot of the function, tested against the
            %   edge tangents (edges only, see jump).
            val = FormValue.jump('Rot','t',pts,pts_ref,parts,c,fe);
        end

        function val = jumpdivd2wn_S(pts,pts_ref,parts,c,fe)
            %jumpdivd2wn_S - Jump of the normal derivative of the Laplacian
            %   VAL = jumpdivd2wn_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   jump of div(D^2 u)*n = grad(laplace(u))*n of the scalar
            %   function across the sides PARTS (see jump).
            val = FormValue.jump('DivD2','n',pts,pts_ref,parts,c,fe);
        end

        function val = jumpdivd2wt_S(pts,pts_ref,parts,c,fe)
            %jumpdivd2wt_S - Jump of the tangential derivative of laplace(u)
            %   VAL = jumpdivd2wt_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   jump of div(D^2 u)*t = grad(laplace(u))*t of the scalar
            %   function across the edges PARTS (edges only, see jump).
            val = FormValue.jump('DivD2','t',pts,pts_ref,parts,c,fe);
        end

        function val = jumpd2wt_S(pts,pts_ref,parts,c,fe)
            %jumpd2wt_S - Tangential jump of the Hessian
            %   VAL = jumpd2wt_S(PTS,PTS_REF,PARTS,C,FE) returns the jump
            %   of the Hessian of the scalar function times the edge
            %   tangents, [D^2 u]*t (edges only, see jump).
            val = FormValue.jump('D2','t',pts,pts_ref,parts,c,fe);
        end

        function val = jumpsymcurlwt_S(pts,pts_ref,parts,c,fe)
            %jumpsymcurlwt_S - Tangential jump of the symmetric curl
            %   VAL = jumpsymcurlwt_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   jump of the symmetric part of the curl of the 2d vector
            %   function (see symcurlw_T) times the edge tangents (edges
            %   only, see jump).
            val = FormValue.jump('SymCurl','t',pts,pts_ref,parts,c,fe);
        end

        function val = jumpdivw_S(pts,pts_ref,parts,c,fe)
            %jumpdivw_S - Jump of the divergence
            %   VAL = jumpdivw_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   difference of the divergence of the function on the two
            %   neighboring elements of each side in PARTS: div(u) for
            %   vector and div(s) for matrix-valued functions (see jump).
            val = FormValue.jump('Div','',pts,pts_ref,parts,c,fe);
        end

        function val = jumpdivwn_S(pts,pts_ref,parts,c,fe)
            %jumpdivwn_S - Normal jump of the divergence
            %   VAL = jumpdivwn_S(PTS,PTS_REF,PARTS,C,FE) returns the jump
            %   of the divergence of the function times the side normals
            %   (see jump).
            val = FormValue.jump('Div','n',pts,pts_ref,parts,c,fe);
        end

        function val = jumpdivwt_S(pts,pts_ref,parts,c,fe)
            %jumpdivwt_S - Tangential jump of the divergence
            %   VAL = jumpdivwt_S(PTS,PTS_REF,PARTS,C,FE) returns the jump
            %   of the divergence of the function times the edge tangents
            %   (edges only, see jump).
            val = FormValue.jump('Div','t',pts,pts_ref,parts,c,fe);
        end

        function val = jumpd3wnnn_S(pts,pts_ref,parts,c,fe)
            %jumpd3wnnn_S - Third normal jump of the function
            %   VAL = jumpd3wnnn_S(PTS,PTS_REF,PARTS,C,FE) returns the jump of
            %   the third derivative of the scalar function, tested three times
            %   against the normals of the first neighboring element, from the
            %   third derivatives stored as [xxx yyy xxy xyy] in 2d and
            %   [xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz] in 3d.
            normals = Form.directions('n',pts_ref,parts,fe.mesh);
            normals = Form.getT1(normals{1});
            val = FormValue.sideValues('D3',pts,pts_ref,parts,c,fe);
            val = FormValue.getJump(val);
            if strcmp(Form.sideOf(fe.mesh),'Edge')
                normals = [normals.^3,3*prod(normals(:,[1 1 2]),2),3*prod(normals(:,[1 2 2]),2)];
            else
                n1 = normals(:,1); n2 = normals(:,2); n3 = normals(:,3);
                normals = [normals.^3,3*n1.^2.*n2,3*n1.^2.*n3,3*n1.*n2.^2,3*n2.^2.*n3, ...
                    3*n1.*n3.^2,3*n2.*n3.^2,6*n1.*n2.*n3];
            end
            val = sum(val.*normals,2);
        end

        function val = jumpgradwntt_S(pts,pts_ref,parts,c,fe)
            %jumpgradwntt_S - Jump of the tangential derivative of n'*w*t
            %   VAL = jumpgradwntt_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   jump of the tangential derivative d_t(n'*w*t) of the
            %   symmetric matrix-valued function w with the edge normals
            %   n and the tangents t = (-n_2,n_1) (edges of 2d meshes
            %   only).
            Form.checkSides(fe.mesh,'planarEdges','FormValue','jumpgradwntt_S');
            val = FormValue.sideValues('Grad',pts,pts_ref,parts,c,fe);
            val = FormValue.getJump(val);
            normals = Form.getNormal4ed(pts_ref,parts,fe.mesh);
            normals = Form.getT1(normals);
            n1 = normals(:,1); n2 = normals(:,2);
            val = val(:,:).*[n1.*n2.^2,-n1.*n2.^2,n2.^3-n1.^2.*n2,-n1.^2.*n2,n1.^2.*n2,-n1.*n2.^2+n1.^3];
            val = sum(val,2);
        end

        function val = jumpwxn_S(pts,pts_ref,parts,c,fe)
            %jumpwxn_S - Jump times the side normals
            %   VAL = jumpwxn_S(PTS,PTS_REF,PARTS,C,FE) returns the jump
            %   of the function times the normals of the first
            %   neighboring element of the sides PARTS: on edges the
            %   products [[u]] n of scalars and the outer products
            %   [[u]] n' of vectors, stored flat with the entry (i,j) as
            %   component i+dim*(j-1), and on faces the cross products
            %   [[u]] x n of 3d vectors.
            normals = Form.directions('n',pts_ref,parts,fe.mesh);
            normals = Form.getT1(normals{1});
            val = FormValue.sideValues('',pts,pts_ref,parts,c,fe);
            val = FormValue.getJump(val);
            if strcmp(Form.sideOf(fe.mesh),'Edge')
                val = val.*permute(normals,[1 3 2]);
                val = val(:,:);
            else
                val = cross(val,normals,2);
            end
        end

        function val = jumpnxcurlwxn_S(pts,pts_ref,parts,c,fe,mu)
            %jumpnxcurlwxn_S - Jump of n x (MU*curl(u)) x n
            %   VAL = jumpnxcurlwxn_S(PTS,PTS_REF,PARTS,C,FE,MU) returns the jump
            %   of the tangential part n x (MU*curl(u)) x n of the 3d vector
            %   function u across the faces PARTS with the face normals n (faces
            %   of 3d meshes only). MU is a function handle of the points that
            %   returns one value per point; it is evaluated at the element
            %   barycenters, i.e. taken as piecewise constant. MU is required
            %   (an empty MU raises an error).
            Form.checkSides(fe.mesh,'faces','FormValue','jumpnxcurlwxn_S');
            val = FormValue.sideValues('Curl',pts,pts_ref,parts,c,fe);
            mu4e = mu(fe.mesh.point4e([1/4,1/4,1/4]));
            mu4f = zeros(size(val,1),2);
            mu4f(:,1) = mu4e(fe.mesh.e4f(parts,1));
            mu4f(fe.mesh.e4f(parts,2)>0,2) = mu4e(fe.mesh.e4f(parts(fe.mesh.e4f(parts,2)>0),2));
            val = mu4f(:,1).*val(:,:,1)-mu4f(:,2).*val(:,:,2);
            normals = fe.mesh.normal4f;
            normals = normals(parts,:);
            val = cross(val,normals,2);
            val = cross(normals,val,2);
        end

        function val = averagecurlw_S(pts,pts_ref,parts,c,fe)
            %averagecurlw_S - Mean of the curl across the sides
            %   VAL = averagecurlw_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   mean of the curl of the function over the two
            %   neighboring elements of each side in PARTS.
            val = FormValue.average('Curl','',pts,pts_ref,parts,c,fe);
        end

        function val = averagewtt_S(pts,pts_ref,parts,c,fe)
            %averagewtt_S - Mean of the tangential-tangential value
            %   VAL = averagewtt_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   mean of the matrix-valued function tested twice against
            %   the edge tangents (edges only).
            val = FormValue.average('','tt',pts,pts_ref,parts,c,fe);
        end

        function val = averagegradwtt_S(pts,pts_ref,parts,c,fe)
            %averagegradwtt_S - Mean of the tangential gradient
            %   VAL = averagegradwtt_S(PTS,PTS_REF,PARTS,C,FE) returns the mean
            %   of t'*Du*t, the gradient matrix of the vector function u tested
            %   twice against the edge tangents t (edges only).
            val = FormValue.average('Grad','tt',pts,pts_ref,parts,c,fe);
        end

        function val = averaged2wt_S(pts,pts_ref,parts,c,fe)
            %averaged2wt_S - Mean of the Hessian times the tangent
            %   VAL = averaged2wt_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   mean of the Hessian of the scalar function times the edge
            %   tangents, {D^2 u}*t (edges only).
            val = FormValue.average('D2','t',pts,pts_ref,parts,c,fe);
        end

        function val = averagedivd2wn_S(pts,pts_ref,parts,c,fe)
            %averagedivd2wn_S - Mean normal derivative of the Laplacian
            %   VAL = averagedivd2wn_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   mean of div(D^2 u)*n = grad(laplace(u))*n of the scalar
            %   function over the two neighboring elements of each side in
            %   PARTS.
            val = FormValue.average('DivD2','n',pts,pts_ref,parts,c,fe);
        end

        function val = averagedivd2wt_S(pts,pts_ref,parts,c,fe)
            %averagedivd2wt_S - Mean tangential derivative of the Laplacian
            %   VAL = averagedivd2wt_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   mean of div(D^2 u)*t = grad(laplace(u))*t of the scalar
            %   function over the two neighboring elements of each edge in
            %   PARTS (edges only).
            val = FormValue.average('DivD2','t',pts,pts_ref,parts,c,fe);
        end

        function val = averagedivwn_S(pts,pts_ref,parts,c,fe)
            %averagedivwn_S - Mean of the divergence times the normal
            %   VAL = averagedivwn_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   mean of the divergence of the function times the side
            %   normals over the two neighboring elements of each side in
            %   PARTS.
            val = FormValue.average('Div','n',pts,pts_ref,parts,c,fe);
        end

        function val = averagedivwt_S(pts,pts_ref,parts,c,fe)
            %averagedivwt_S - Mean of the divergence times the tangent
            %   VAL = averagedivwt_S(PTS,PTS_REF,PARTS,C,FE) returns the
            %   mean of the divergence of the function times the edge
            %   tangents over the two neighboring elements of each edge in
            %   PARTS (edges only).
            val = FormValue.average('Div','t',pts,pts_ref,parts,c,fe);
        end

        %% side values
        function val = sideValues(name,pts,pts_ref,parts,c,fe)
            %sideValues - Function values on both neighboring elements
            %   VAL = sideValues(NAME,PTS,PTS_REF,PARTS,C,FE) returns the
            %   values of the finite element function with coefficients
            %   C of FE on the two neighboring elements of the sides
            %   PARTS (see Form.sideOf), with the elements along the
            %   last dimension; the values of the missing second element of a
            %   boundary side are zero. NAME selects the derivative: ''
            %   (values), 'Grad', 'D2', 'D3', 'Rot', 'Curl' or 'Div'.
            %   FE provides them directly by getEdge<NAME>Value or
            %   getFace<NAME>Value, or they are reconstructed from its
            %   element values get<NAME>Value. Three values are derived
            %   from these: 'SymGrad', the symmetric part of the gradient
            %   with the components [11 22 12] or [11 22 33 12 13 23],
            %   'DivD2', div(D^2 u) = grad(laplace(u)) of scalar functions
            %   from 'D3' (see divd2w_T), and 'SymCurl', the symmetric part
            %   of the row-wise curl of 2d vector functions from 'Grad'
            %   (see symcurlw_T).
            switch name
                case 'SymGrad'
                    val = FormValue.sideValues('Grad',pts,pts_ref,parts,c,fe);
                    val = FormValue.typeOf(val,fe).symPart(val,2);
                    return
                case 'DivD2'
                    d3 = FormValue.sideValues('D3',pts,pts_ref,parts,c,fe);
                    if size(d3,2) == 4
                        val = [sum(d3(:,[1 4],:),2),sum(d3(:,[2 3],:),2)];
                    else
                        val = [sum(d3(:,[1 6 8],:),2),sum(d3(:,[4 2 9],:),2),sum(d3(:,[5 7 3],:),2)];
                    end
                    return
                case 'SymCurl'
                    grad = FormValue.sideValues('Grad',pts,pts_ref,parts,c,fe);
                    if size(fe.mesh.c4n,2) ~= 2 || ndims(grad) ~= 3 || size(grad,2) ~= 4
                        error('fem4m:FormValue:symCurlNotDefined', ...
                            'FormValue: sym(curl(u)) is only defined for 2d vector fields');
                    end
                    % grad = [dx v1, dx v2, dy v1, dy v2], curl column by column
                    curl = [grad(:,3:4,:),-grad(:,1:2,:)];
                    val = FormType.matrix2d.symPart(curl,2);
                    return
            end
            side = Form.sideOf(fe.mesh);
            if strcmp(side,'Edge')
                elemValues = @getEdgeValues4ElemFunktion;
                refValues = @getEdgeValues4RefElemFunktion;
            else
                elemValues = @getFaceValues4ElemFunktion;
                refValues = @getFaceValues4RefElemFunktion;
            end
            sideMethod = ['get',side,name,'Value'];
            elemMethod = ['get',name,'Value'];
            if ismethod(fe,sideMethod)
                val = fe.(sideMethod)(pts,pts_ref,parts,c);
            elseif ~fe.isRefElemBased
                val = elemValues(pts,parts,@(x,ind)fe.(elemMethod)(x,[],ind,c),fe.mesh);
            else
                val = refValues(pts_ref,@(x)fe.(elemMethod)([],x,[],c),fe.mesh);
                val = val(parts,:,:,:,:);
            end
        end

        %% side operators
        function val = oneSide(name,letters,pts,pts_ref,parts,c,fe,alpha)
            %oneSide - Side operator on the first neighboring element
            %   VAL = oneSide(NAME,LETTERS,PTS,PTS_REF,PARTS,C,FE)
            %   evaluates the values NAME (see sideValues) of the finite
            %   element function with coefficients C of FE on the first
            %   neighboring element of the sides PARTS and multiplies
            %   them by the directions LETTERS of that element: none for
            %   '', the normal 'n' or tangent 't' (see oneSideTimes), or
            %   for two letters AB the matrix values between them,
            %   B'*M*A (see contractOneSide).
            %
            %   VAL = oneSide(...,ALPHA) applies the coefficient ALPHA to the
            %   values before the directions: for 'Grad' a matrix acting on all
            %   gradient components, for 'SymGrad' the elasticity tensor acting
            %   on the symmetric gradient and for the other NAMEs a scalar for
            %   scalar values and a d-by-d matrix acting on vectors and on the
            %   columns of matrices (see FormEval.evalCoefficient).
            dirs = Form.directions(letters,pts_ref,parts,fe.mesh);
            dirs = cellfun(@Form.getT1,dirs,'UniformOutput',false);
            val = Form.getT1(FormValue.sideValues(name,pts,pts_ref,parts,c,fe));
            if nargin > 7
                val = FormValue.applyCoefficient(name,pts,pts_ref,parts,val,fe,alpha);
            end
            val = FormValue.timesDirections(val,dirs,fe);
        end

        function val = jump(name,letters,pts,pts_ref,parts,c,fe,alpha)
            %JUMP - Jump of a side operator across the sides
            %   VAL = JUMP(NAME,LETTERS,PTS,PTS_REF,PARTS,C,FE) returns
            %   the jump of the values NAME (see sideValues) of the
            %   finite element function with coefficients C of FE times
            %   the directions LETTERS (see oneSide) across the sides
            %   PARTS. Each element is multiplied by its own directions,
            %   which on the second element are the negated ones, or
            %   its conormals on surfaces. With one direction the sum of
            %   both elements is the jump; with none or two the signs
            %   cancel and the second element is subtracted. On a boundary
            %   side the jump is the value of the first element.
            %
            %   VAL = JUMP(...,ALPHA) applies the coefficient ALPHA to the
            %   values of both elements before the directions (see
            %   oneSide).
            dirs = Form.directions(letters,pts_ref,parts,fe.mesh);
            val = FormValue.sideValues(name,pts,pts_ref,parts,c,fe);
            sides = cell(1,2);
            for s = 1 : 2
                curVal = FormValue.sideSlice(val,s);
                if nargin > 7
                    curVal = FormValue.applyCoefficient(name,pts,pts_ref,parts,curVal,fe,alpha);
                end
                curDirs = cellfun(@(d)d(:,:,s),dirs,'UniformOutput',false);
                sides{s} = FormValue.timesDirections(curVal,curDirs,fe);
            end
            if mod(numel(dirs),2) == 0
                val = sides{1} - sides{2};
            else
                val = sides{1} + sides{2};
            end
        end

        function val = average(name,letters,pts,pts_ref,parts,c,fe,alpha)
            %AVERAGE - Mean of a side operator across the sides
            %   VAL = AVERAGE(NAME,LETTERS,PTS,PTS_REF,PARTS,C,FE)
            %   returns the mean of the values NAME (see sideValues) of
            %   the finite element function with coefficients C of FE
            %   over the two neighboring elements of the sides PARTS
            %   (see getAverage) times the directions LETTERS of the
            %   first element (see oneSide). On a boundary side the mean is
            %   the value of the first element.
            %
            %   VAL = AVERAGE(...,ALPHA) applies the coefficient ALPHA to
            %   the mean before the directions (see oneSide).
            e4s = Form.sideOf(fe.mesh,'e4s');
            dirs = Form.directions(letters,pts_ref,parts,fe.mesh);
            dirs = cellfun(@Form.getT1,dirs,'UniformOutput',false);
            val = FormValue.sideValues(name,pts,pts_ref,parts,c,fe);
            val = FormValue.getAverage(val,e4s(parts,:));
            if nargin > 7
                val = FormValue.applyCoefficient(name,pts,pts_ref,parts,val,fe,alpha);
            end
            val = FormValue.timesDirections(val,dirs,fe);
        end

        %% value type
        function type = typeOf(val,fe)
            %typeOf - Value type of function values
            %   TYPE = typeOf(VAL,FE) returns the FormType of the values
            %   VAL of a finite element function of FE, stored as rows
            %   x components (x 2 neighboring elements for side values).
            type = FormType.of(size(val,2),size(fe.mesh.c4n,2));
        end

        %% products with directions
        function val = oneSideTimes(val,dirs,fe)
            %oneSideTimes - One-sided values times directions
            %   VAL = oneSideTimes(VAL,DIRS,FE) multiplies the values
            %   VAL of FE (rows x components, see typeOf) by the
            %   directions DIRS (rows x dim): scalars times DIRS, the
            %   inner product for vectors and the matrix times DIRS
            %   for matrices.
            type = FormValue.typeOf(val,fe);
            if type.isMatrix()
                val = sum(type.toMatrix(val,2).*permute(dirs,[1 3 2]),3);
            elseif type.isVector()
                val = sum(val.*dirs,2);
            else
                val = val.*dirs;
            end
        end

        function val = contractOneSide(val,left,right,fe)
            %contractOneSide - One-sided matrix values between two directions
            %   VAL = contractOneSide(VAL,LEFT,RIGHT,FE) returns
            %   LEFT'*M*RIGHT for the matrix values M in VAL of FE
            %   (rows x components, see typeOf) and the directions
            %   LEFT and RIGHT (rows x dim).
            weights = left.*permute(right,[1 3 2]);
            val = sum(FormValue.typeOf(val,fe).toMatrix(val,2).*weights,[2 3]);
        end

        %% jumps & averages (value)
        function val = getJump(val)
            %getJump - Difference of the two neighboring values
            %   VAL = getJump(VAL) returns the values of the first neighboring
            %   element in VAL minus those of the second one, stored along the
            %   last dimension of VAL (see sideValues).
            nDim = numel(size(val));
            val = shiftdim(val,nDim-1);
            originDim = size(val);
            val = val(:,:);
            val = val(1,:)-val(2,:);
            val = reshape(val,[1,originDim(2:end)]);
            val = shiftdim(val,1);
        end

        function val = getAverage(val,e4s)
            %getAverage - Mean of the two neighboring values
            %   VAL = getAverage(VAL,E4S) returns the mean of the values of the
            %   two neighboring elements in VAL, stored along the last
            %   dimension (see sideValues). E4S holds the two neighboring
            %   elements of the side of each row of VAL; where its second
            %   entry is zero (boundary side), the mean is doubled to the
            %   value of the first element.
            nDim = numel(size(val));
            val = shiftdim(val,nDim-1);
            originDim = size(val);
            val = val(:,:);
            val = 0.5*(val(1,:)+val(2,:));
            val = reshape(val,[1,originDim(2:end)]);
            val = shiftdim(val,1);
            % value on the boundary
            originDim = size(val);
            val(e4s(:,2)==0,:) = 2.*val(e4s(:,2)==0,:);
            val = reshape(val,originDim);
        end

    end

    methods (Static, Access = private)

        function val = applyCoefficient(name,pts,pts_ref,parts,val,fe,alpha)
            %applyCoefficient - Coefficient applied to one-sided values
            %   VAL = applyCoefficient(NAME,PTS,PTS_REF,PARTS,VAL,FE,ALPHA)
            %   applies the coefficient ALPHA to the values NAME of FE on
            %   one side (rows x components, see oneSide for the rules).
            %   ALPHA turns symmetric matrix values into full ones.
            dim = size(fe.mesh.c4n,2);
            type = FormValue.typeOf(val,fe);
            switch name
                case {'Grad','SymGrad'}
                    alpha = FormEval.evalCoefficient(pts,pts_ref,parts,alpha,type);
                    val = matMul(alpha,val);
                otherwise
                    if type.isScalar()
                        alpha = FormEval.evalCoefficient(pts,pts_ref,parts,alpha,FormType.scalar);
                        val = alpha.*val;
                        return
                    end
                    alpha = FormEval.evalCoefficient(pts,pts_ref,parts,alpha,FormType.vector(dim));
                    if type.isMatrix()
                        val = matMul(alpha,type.toMatrix(val,2));
                        val = val(:,:);
                    else
                        val = matMul(alpha,val);
                    end
            end
        end

        function val = timesDirections(val,dirs,fe)
            %timesDirections - One-sided values times their directions
            %   VAL = timesDirections(VAL,DIRS,FE) keeps the one-sided
            %   values VAL of FE without directions, multiplies them by
            %   one direction (see oneSideTimes) or returns the matrix
            %   values between two directions, DIRS{2}'*M*DIRS{1} (see
            %   contractOneSide).
            if isscalar(dirs)
                val = FormValue.oneSideTimes(val,dirs{1},fe);
            elseif numel(dirs) == 2
                val = FormValue.contractOneSide(val,dirs{2},dirs{1},fe);
            end
        end

        function val = sideSlice(val,s)
            %sideSlice - Values of one side
            %   VAL = sideSlice(VAL,S) returns the values of the side S of
            %   the two-sided values VAL, stored along the last dimension.
            idx = repmat({':'},1,ndims(val));
            idx{end} = s;
            val = val(idx{:});
        end

    end
end
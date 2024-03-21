classdef FormBasis
    %FormBasis - Basis function operators for form assembly
    %   FormBasis provides static methods that evaluate differential
    %   operators of finite element basis functions on the elements
    %   (suffix _T) and on the sides (suffix _S), i.e. the edges of 2d,
    %   polygonal and surface meshes and the faces of 3d meshes (see
    %   Form.sideOf), including jump and average variants. BilinearForm
    %   and LinearForm select them by name from the parsed form string.
    %
    %   Every operator takes the points PTS (one row per element or side of
    %   PARTS), the reference point PTS_REF of the element or side (1-by-d, d
    %   the dimension of the reference element or side, so 1-by-1 on edges;
    %   empty for polygonal elements), the indices PARTS of the elements
    %   (suffix _T) or sides (suffix _S) and the finite element FE, as
    %   integrate passes them through the forms (see FE for which of PTS and
    %   PTS_REF an element uses). VAL has one row per element or side of PARTS,
    %   the basis functions along dim 2 and the components from dim 3 on (see
    %   typeOf). The jump and average operators return the basis functions of
    %   both neighboring elements, those of the second following those of the
    %   first along dim 2 (zero for the missing neighbor of a boundary side);
    %   the other side operators use the first neighboring element only. The
    %   tangential operators (letter t) are defined on edges only; on faces
    %   they raise the error fem4m:Form:noFaceTangent.
    %
    %   FormBasis functions:
    %       w_T              - Basis functions on the element
    %       wn_T             - Normal component of the basis
    %       trw_T            - Trace of the basis
    %       trgradw_T        - Trace of the basis gradient
    %       gradw_T          - Gradient of the basis
    %       curlw_T          - Curl of the basis
    %       rotw_T           - Rot of the basis
    %       divw_T           - Divergence of the basis
    %       epsw_T           - Symmetric gradient of the basis
    %       devw_T           - Deviator of the basis
    %       symw_T           - Symmetric part of the basis
    %       d2w_T            - Hessian of the basis
    %       laplacew_T       - Laplacian of the basis
    %       curlcurlw_T      - Curl-curl of the basis
    %       curlrotw_T       - Curl of the rot of the basis
    %       rotrotw_T        - Rot of the rot of the basis
    %       symcurlw_T       - Symmetric part of the curl of the basis
    %       divdivw_T        - Double divergence of the basis
    %       divd2w_T         - Divergence of the Hessian of the basis
    %       w_S              - Basis functions on the sides
    %       gradw_S          - Gradient of the basis on the sides
    %       d2w_S            - Hessian of the basis on the sides
    %       curlw_S          - Curl of the basis on the sides
    %       wn_S             - Normal component of the basis on the sides
    %       wt_S             - Tangential component of the basis on the edges
    %       wnn_S            - Normal-normal component of the basis
    %       d2wnn_S          - Second normal derivative on the sides
    %       d2wn_S           - Hessian times the normal on the sides
    %       gradwn_S         - Normal derivative of the basis on the sides
    %       divwn_S          - Divergence times the normal on the sides
    %       divd2wn_S        - Normal derivative of the Laplacian on the sides
    %       gradwxn_S        - Gradient crossed with the face normal
    %       jumpw_S          - Jump of the basis across the sides
    %       jumpgradw_S      - Jump of the basis gradient
    %       jumpwn_S         - Normal jump of the basis across the sides
    %       jumpwt_S         - Tangential jump of the basis across the edges
    %       jumpgradwn_S     - Normal jump of the basis gradient
    %       jumpwnn_S        - Jump of the normal-normal component
    %       jumpwtn_S        - Jump of the tangential-normal component
    %       jumpwnt_S        - Jump of the normal-tangential component
    %       jumpwxn_S        - Jump times the side normals
    %       jumpgradwxn_S    - Cross jump of the basis gradient
    %       jumpsymcurlwt_S  - Tangential jump of the symmetric curl
    %       jumpd2wt_S       - Tangential jump of the Hessian
    %       jumpdivd2wt_S    - Tangential jump of the gradient of the Laplacian
    %       jumpdivw_S       - Jump of the divergence
    %       jumpdivwn_S      - Normal jump of the divergence
    %       jumpdivwt_S      - Tangential jump of the divergence
    %       averagew_S       - Mean of the basis across the sides
    %       averagegradw_S   - Mean of the basis gradient
    %       averagecurlw_S   - Mean of the curl across the sides
    %       averagewn_S      - Mean of the normal basis
    %       averagewt_S      - Mean of the tangential basis
    %       averagegradwn_S  - Mean of the normal basis gradient
    %       averagegradwt_S  - Mean of the tangential basis gradient
    %       averagewnn_S     - Mean of the normal-normal component
    %       averaged2wn_S    - Mean of the Hessian times the normal
    %       averaged2wt_S    - Mean of the Hessian times the tangent
    %       averaged2wnn_S   - Mean of the second normal derivative
    %       averagedivd2wn_S - Mean of the normal derivative of the Laplacian
    %       averagedivd2wt_S - Mean tangential derivative of the Laplacian
    %       averagedivwn_S   - Mean of the divergence times the normal
    %       averagedivwt_S   - Mean of the divergence times the tangent
    %       sideValues       - Basis values on both neighboring elements
    %       oneSide          - Side operator on the first neighboring element
    %       jump             - Jump of a side operator across the sides
    %       average          - Mean of a side operator across the sides
    %       typeOf           - Value type of element basis values
    %       sideTypeOf       - Value type of two-sided basis values
    %       oneSideTimes     - One-sided basis values times directions
    %       contractOneSide  - One-sided matrix values between two directions
    %       sideTimes        - Two-sided basis values times directions
    %       contractSides    - Two-sided matrix values between two directions
    %       elementRows      - Keep the rows of the elements PARTS
    %       getJump          - Negate the second side values
    %       stackSides       - Stack the two side values along the basis
    %       getAverage       - Average side values of the basis
    %
    %   See also FormValue, Form, BilinearForm, LinearForm, FormType

    % Copyright (c) 2026 by the fem4m authors

    methods(Static)

        %% basis on volume
        function val = w_T(pts,pts_ref,parts,fe)
            %w_T - Basis functions on the element
            %   VAL = w_T(PTS,PTS_REF,PARTS,FE) evaluates the basis
            %   functions of FE in the elements PARTS at the reference
            %   point PTS_REF (1-by-d), or at the points PTS (one row per
            %   element of PARTS) for elements that are not
            %   reference-element based (see FE).
            val = FormBasis.elementRows(fe.getBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = wn_T(pts,pts_ref,parts,fe)
            %wn_T - Normal component of the basis
            %   VAL = wn_T(PTS,PTS_REF,PARTS,FE) evaluates the normal
            %   component of the basis functions of FE, for elements
            %   that provide getNormalBasis (no element does yet).
            val = FormBasis.elementRows(fe.getNormalBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = trw_T(pts,pts_ref,parts,fe)
            %trw_T - Trace of the basis
            %   VAL = trw_T(PTS,PTS_REF,PARTS,FE) evaluates the trace
            %   of the matrix-valued basis of FE by summing its diagonal
            %   entries, or by its getTraceBasis if FE provides one (no
            %   element does yet).
            if ismethod(fe,'getTraceBasis')
                val = FormBasis.elementRows(fe.getTraceBasis(pts,pts_ref,parts),parts,fe);
            else
                val = FormBasis.elementRows(fe.getBasis(pts,pts_ref,parts),parts,fe);
                type = FormBasis.typeOf(val,fe);
                val = sum(val(:,:,type.diagonal),3);
            end
        end

        function val = trgradw_T(pts,pts_ref,parts,fe)
            %trgradw_T - Trace of the basis gradient
            %   VAL = trgradw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   trace of the gradient of the basis of FE, i.e. the
            %   divergence of a vector basis.
            val = FormBasis.elementRows(fe.getGradBasis(pts,pts_ref,parts),parts,fe);
            type = FormBasis.typeOf(val,fe);
            val = sum(val(:,:,type.diagonal),3);
        end

        function val = gradw_T(pts,pts_ref,parts,fe)
            %gradw_T - Gradient of the basis
            %   VAL = gradw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   gradient of the basis functions of FE.
            val = FormBasis.elementRows(fe.getGradBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = curlw_T(pts,pts_ref,parts,fe)
            %curlw_T - Curl of the basis
            %   VAL = curlw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   curl of the basis functions of FE: the vector curl
            %   (d_y u, -d_x u) of a 2d scalar basis and the curl of a
            %   3d vector basis.
            val = FormBasis.elementRows(fe.getCurlBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = rotw_T(pts,pts_ref,parts,fe)
            %rotw_T - Rot of the basis
            %   VAL = rotw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   rotation d_x u_2 - d_y u_1 of the 2d vector basis
            %   functions of FE, row-wise for matrix-valued ones.
            val = FormBasis.elementRows(fe.getRotBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = divw_T(pts,pts_ref,parts,fe)
            %divw_T - Divergence of the basis
            %   VAL = divw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   divergence of the basis functions of FE.
            val = FormBasis.elementRows(fe.getDivBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = epsw_T(pts,pts_ref,parts,fe)
            %epsw_T - Symmetric gradient of the basis
            %   VAL = epsw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   symmetric gradient of the basis functions of FE.
            val = FormBasis.elementRows(fe.getSymGradBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = devw_T(pts,pts_ref,parts,fe)
            %devw_T - Deviator of the basis
            %   VAL = devw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   deviator of the basis functions of FE, for elements that
            %   provide getDevBasis (no element does yet).
            val = FormBasis.elementRows(fe.getDevBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = symw_T(pts,pts_ref,parts,fe)
            %symw_T - Symmetric part of the basis
            %   VAL = symw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   symmetric part (M+M')/2 of the matrix-valued basis
            %   functions M of FE, stored as [11 22 12] in 2d and
            %   [11 22 33 12 13 23] in 3d (see FormType.symPart).
            val = FormBasis.elementRows(fe.getBasis(pts,pts_ref,parts),parts,fe);
            val = FormBasis.typeOf(val,fe).symPart(val,3);
        end

        function val = d2w_T(pts,pts_ref,parts,fe)
            %d2w_T - Hessian of the basis
            %   VAL = d2w_T(PTS,PTS_REF,PARTS,FE) evaluates the second
            %   derivatives of the basis functions of FE.
            val = FormBasis.elementRows(fe.getD2Basis(pts,pts_ref,parts),parts,fe);
        end

        function val = laplacew_T(pts,pts_ref,parts,fe)
            %laplacew_T - Laplacian of the basis
            %   VAL = laplacew_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   Laplacian of the basis functions of FE by summing the
            %   diagonal of the Hessian, component by component for
            %   vector- and matrix-valued bases, whose Hessians are
            %   stored along dim 4.
            val = FormBasis.elementRows(fe.getD2Basis(pts,pts_ref,parts),parts,fe);
            if ndims(val) == 4
                val = sum(val(:,:,:,1:size(pts,2)),4);
            else
                val = sum(val(:,:,1:size(pts,2)),3);
            end
        end

        function val = curlcurlw_T(pts,pts_ref,parts,fe)
            %curlcurlw_T - Curl-curl of the basis
            %   VAL = curlcurlw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   curl-curl of the basis functions of FE: the Airy
            %   matrices curl(curl(u)) = [u_yy -u_xy; -u_xy u_xx] of a
            %   2d scalar basis and the curl-curls of a 3d vector basis.
            val = FormBasis.elementRows(fe.getCurlCurlBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = curlrotw_T(pts,pts_ref,parts,fe)
            %curlrotw_T - Curl of the rot of the basis
            %   VAL = curlrotw_T(PTS,PTS_REF,PARTS,FE) evaluates
            %   curl(rot(u)) of the 2d vector basis functions of FE.
            val = FormBasis.elementRows(fe.getCurlRotBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = rotrotw_T(pts,pts_ref,parts,fe)
            %rotrotw_T - Rot of the rot of the basis
            %   VAL = rotrotw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   scalar rot(rot(sigma)) of the 2d symmetric matrix basis
            %   functions of FE.
            val = FormBasis.elementRows(fe.getRotRotBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = symcurlw_T(pts,pts_ref,parts,fe)
            %symcurlw_T - Symmetric part of the curl of the basis
            %   VAL = symcurlw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   symmetric part of the row-wise curl, with the rows
            %   (d_y v_i, -d_x v_i), of the 2d vector basis functions v
            %   of FE: [d_y v_1, -d_x v_2, (d_y v_2 - d_x v_1)/2]. Other
            %   bases raise the error fem4m:FormBasis:symCurlNotDefined.
            grad = FormBasis.elementRows(fe.getGradBasis(pts,pts_ref,parts),parts,fe);
            if size(fe.mesh.c4n,2) ~= 2 || ndims(grad) ~= 3 || size(grad,3) ~= 4
                error('fem4m:FormBasis:symCurlNotDefined', ...
                    'FormBasis: sym(curl(u)) is only defined for 2d vector fields');
            end
            % grad = [dx v1, dx v2, dy v1, dy v2], curl column by column
            curl = cat(3,grad(:,:,3:4),-grad(:,:,1:2));
            val = FormType.matrix2d.symPart(curl,3);
        end

        function val = divdivw_T(pts,pts_ref,parts,fe)
            %divdivw_T - Double divergence of the basis
            %   VAL = divdivw_T(PTS,PTS_REF,PARTS,FE) evaluates the
            %   double divergence sum_ij d_i d_j s_ij of the symmetric
            %   matrix basis functions s of FE.
            val = FormBasis.elementRows(fe.getDivDivBasis(pts,pts_ref,parts),parts,fe);
        end

        function val = divd2w_T(pts,pts_ref,parts,fe)
            %divd2w_T - Divergence of the Hessian of the basis
            %   VAL = divd2w_T(PTS,PTS_REF,PARTS,FE) evaluates div(D^2 v)
            %   = grad(laplace(v)) of the scalar basis functions v of FE
            %   from the third derivatives stored as [xxx yyy xxy xyy] in
            %   2d and [xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz] in 3d.
            d3 = FormBasis.elementRows(fe.getD3Basis(pts,pts_ref,parts),parts,fe);
            if size(d3,3) == 4
                val = cat(3,sum(d3(:,:,[1 4]),3),sum(d3(:,:,[2 3]),3));
            else
                val = cat(3,sum(d3(:,:,[1 6 8]),3),sum(d3(:,:,[4 2 9]),3),sum(d3(:,:,[5 7 3]),3));
            end
        end

        %% basis on sides (edges, or faces of 3d meshes; see Form.sideOf)
        function val = w_S(pts,pts_ref,parts,fe)
            %w_S - Basis functions on the sides
            %   VAL = w_S(PTS,PTS_REF,PARTS,FE) evaluates the basis
            %   functions of FE on the first neighboring element of the
            %   sides PARTS.
            val = FormBasis.oneSide('','',pts,pts_ref,parts,fe);
        end

        function val = gradw_S(pts,pts_ref,parts,fe)
            %gradw_S - Gradient of the basis on the sides
            %   VAL = gradw_S(PTS,PTS_REF,PARTS,FE) evaluates the
            %   gradient of the basis functions of FE on the first
            %   neighboring element of the sides PARTS.
            val = FormBasis.oneSide('Grad','',pts,pts_ref,parts,fe);
        end

        function val = d2w_S(pts,pts_ref,parts,fe)
            %d2w_S - Hessian of the basis on the sides
            %   VAL = d2w_S(PTS,PTS_REF,PARTS,FE) evaluates the second
            %   derivatives of the basis functions of FE on the first
            %   neighboring element of the sides PARTS.
            val = FormBasis.oneSide('D2','',pts,pts_ref,parts,fe);
        end

        function val = curlw_S(pts,pts_ref,parts,fe)
            %curlw_S - Curl of the basis on the sides
            %   VAL = curlw_S(PTS,PTS_REF,PARTS,FE) evaluates the curl
            %   of the basis functions of FE (see curlw_T) on the first
            %   neighboring element of the sides PARTS.
            val = FormBasis.oneSide('Curl','',pts,pts_ref,parts,fe);
        end

        function val = wn_S(pts,pts_ref,parts,fe)
            %wn_S - Normal component of the basis on the sides
            %   VAL = wn_S(PTS,PTS_REF,PARTS,FE) multiplies the basis of
            %   FE on the first neighboring element of the sides PARTS
            %   by the normal of that element.
            val = FormBasis.oneSide('','n',pts,pts_ref,parts,fe);
        end

        function val = wt_S(pts,pts_ref,parts,fe)
            %wt_S - Tangential component of the basis on the edges
            %   VAL = wt_S(PTS,PTS_REF,PARTS,FE) multiplies the basis of
            %   FE on the first neighboring element of the edges PARTS
            %   by the tangent of that element (edges only).
            val = FormBasis.oneSide('','t',pts,pts_ref,parts,fe);
        end

        function val = wnn_S(pts,pts_ref,parts,fe)
            %wnn_S - Normal-normal component of the basis
            %   VAL = wnn_S(PTS,PTS_REF,PARTS,FE) contracts the
            %   matrix-valued basis of FE on the first neighboring
            %   element of the sides PARTS with the outer products of
            %   the normal of that element.
            val = FormBasis.oneSide('','nn',pts,pts_ref,parts,fe);
        end

        function val = d2wnn_S(pts,pts_ref,parts,fe)
            %d2wnn_S - Second normal derivative on the sides
            %   VAL = d2wnn_S(PTS,PTS_REF,PARTS,FE) contracts the second
            %   derivatives of the basis of FE on the first neighboring
            %   element of the sides PARTS with the outer products of
            %   the normal of that element.
            val = FormBasis.oneSide('D2','nn',pts,pts_ref,parts,fe);
        end

        function val = d2wn_S(pts,pts_ref,parts,fe)
            %d2wn_S - Hessian times the normal on the sides
            %   VAL = d2wn_S(PTS,PTS_REF,PARTS,FE) multiplies the second
            %   derivatives of the scalar basis of FE on the first
            %   neighboring element of the sides PARTS by the normal of
            %   that element, D^2 v*n.
            val = FormBasis.oneSide('D2','n',pts,pts_ref,parts,fe);
        end

        function val = gradwn_S(pts,pts_ref,parts,fe)
            %gradwn_S - Normal derivative of the basis on the sides
            %   VAL = gradwn_S(PTS,PTS_REF,PARTS,FE) multiplies the
            %   gradient of the basis of FE on the first neighboring
            %   element of the sides PARTS by the normal of that
            %   element.
            val = FormBasis.oneSide('Grad','n',pts,pts_ref,parts,fe);
        end

        function val = divwn_S(pts,pts_ref,parts,fe)
            %divwn_S - Divergence times the normal on the sides
            %   VAL = divwn_S(PTS,PTS_REF,PARTS,FE) multiplies the
            %   divergence of the basis of FE on the first neighboring
            %   element of the sides PARTS by the normal of that element:
            %   div(v)*n for vector and div(s)*n for matrix-valued bases.
            val = FormBasis.oneSide('Div','n',pts,pts_ref,parts,fe);
        end

        function val = divd2wn_S(pts,pts_ref,parts,fe)
            %divd2wn_S - Normal derivative of the Laplacian on the sides
            %   VAL = divd2wn_S(PTS,PTS_REF,PARTS,FE) returns div(D^2 v)*n
            %   = grad(laplace(v))*n of the scalar basis of FE on the first
            %   neighboring element of the sides PARTS with the normal of
            %   that element.
            val = FormBasis.oneSide('DivD2','n',pts,pts_ref,parts,fe);
        end

        function val = gradwxn_S(pts,pts_ref,parts,fe)
            %gradwxn_S - Gradient crossed with the face normal
            %   VAL = gradwxn_S(PTS,PTS_REF,PARTS,FE) returns the cross
            %   product of the gradient of the basis of FE on the first
            %   neighboring element of the faces PARTS with the face
            %   normals (faces of 3d meshes only; other sides raise the
            %   error fem4m:FormBasis:sideDimension).
            Form.checkSides(fe.mesh,'faces','FormBasis','gradwxn_S');
            val = FormBasis.sideValues('Grad',pts,pts_ref,parts,fe);
            val = Form.getT1(val);
            normals = fe.mesh.normal4f;
            normals = permute(normals(parts,:),[1 3 2]).*ones(1,size(val,2));
            val = cross(val,normals,3);
        end

        %% jumps & averages (basis)
        function val = jumpw_S(pts,pts_ref,parts,fe)
            %jumpw_S - Jump of the basis across the sides
            %   VAL = jumpw_S(PTS,PTS_REF,PARTS,FE) returns the
            %   difference of the basis values of FE on the two
            %   neighboring elements of the sides PARTS, stacked along
            %   dim 2.
            val = FormBasis.jump('','',pts,pts_ref,parts,fe);
        end

        function val = jumpgradw_S(pts,pts_ref,parts,fe)
            %jumpgradw_S - Jump of the basis gradient
            %   VAL = jumpgradw_S(PTS,PTS_REF,PARTS,FE) returns the
            %   difference of the gradient of the basis of FE on the two
            %   neighboring elements of the sides PARTS, stacked along
            %   dim 2.
            val = FormBasis.jump('Grad','',pts,pts_ref,parts,fe);
        end

        function val = jumpwn_S(pts,pts_ref,parts,fe)
            %jumpwn_S - Normal jump of the basis across the sides
            %   VAL = jumpwn_S(PTS,PTS_REF,PARTS,FE) returns the jump of
            %   the basis of FE times the normals across the sides
            %   PARTS (see jump).
            val = FormBasis.jump('','n',pts,pts_ref,parts,fe);
        end

        function val = jumpwt_S(pts,pts_ref,parts,fe)
            %jumpwt_S - Tangential jump of the basis across the edges
            %   VAL = jumpwt_S(PTS,PTS_REF,PARTS,FE) returns the jump of
            %   the basis of FE times the tangents across the edges
            %   PARTS (edges only, see jump).
            val = FormBasis.jump('','t',pts,pts_ref,parts,fe);
        end

        function val = jumpgradwn_S(pts,pts_ref,parts,fe)
            %jumpgradwn_S - Normal jump of the basis gradient
            %   VAL = jumpgradwn_S(PTS,PTS_REF,PARTS,FE) returns the jump
            %   of the normal derivative of the basis of FE across the
            %   sides PARTS (see jump).
            val = FormBasis.jump('Grad','n',pts,pts_ref,parts,fe);
        end

        function val = jumpwnn_S(pts,pts_ref,parts,fe)
            %jumpwnn_S - Jump of the normal-normal component
            %   VAL = jumpwnn_S(PTS,PTS_REF,PARTS,FE) returns the jump of
            %   the normal-normal components of the matrix-valued basis
            %   of FE across the sides PARTS (see jump).
            val = FormBasis.jump('','nn',pts,pts_ref,parts,fe);
        end

        function val = jumpwtn_S(pts,pts_ref,parts,fe)
            %jumpwtn_S - Jump of the tangential-normal component
            %   VAL = jumpwtn_S(PTS,PTS_REF,PARTS,FE) returns the jump of
            %   n'*M*t of the matrix-valued basis M of FE across the
            %   edges PARTS (edges only, see jump).
            val = FormBasis.jump('','tn',pts,pts_ref,parts,fe);
        end

        function val = jumpwnt_S(pts,pts_ref,parts,fe)
            %jumpwnt_S - Jump of the normal-tangential component
            %   VAL = jumpwnt_S(PTS,PTS_REF,PARTS,FE) returns the jump of
            %   t'*M*n of the matrix-valued basis M of FE across the
            %   edges PARTS (edges only, see jump).
            val = FormBasis.jump('','nt',pts,pts_ref,parts,fe);
        end

        function val = jumpwxn_S(pts,pts_ref,parts,fe)
            %jumpwxn_S - Jump times the side normals
            %   VAL = jumpwxn_S(PTS,PTS_REF,PARTS,FE) returns the jump of
            %   the basis of FE times the normals of the first
            %   neighboring element of the sides PARTS: on edges the
            %   products [[v]] n of scalars and the outer products
            %   [[v]] n' of vectors, stored flat with the entry (i,j) as
            %   component i+dim*(j-1), and on faces the cross products
            %   [[v]] x n of 3d vectors.
            normals = Form.directions('n',pts_ref,parts,fe.mesh);
            normals = Form.getT1(normals{1});
            val = FormBasis.sideValues('',pts,pts_ref,parts,fe);
            val = FormBasis.getJump(val);
            val = FormBasis.stackSides(val);
            if strcmp(Form.sideOf(fe.mesh),'Edge')
                val = permute(normals,[1 3 4 2]).*val;
                val = val(:,:,:);
            else
                normals = permute(normals,[1 3 2]).*ones(1,size(val,2));
                val = cross(val,normals,3);
            end
        end

        function val = jumpgradwxn_S(pts,pts_ref,parts,fe)
            %jumpgradwxn_S - Cross jump of the basis gradient
            %   VAL = jumpgradwxn_S(PTS,PTS_REF,PARTS,FE) returns the
            %   cross products of the jump of the gradient of the basis
            %   of FE with the face normals (faces of 3d meshes only;
            %   other sides raise the error fem4m:FormBasis:sideDimension).
            Form.checkSides(fe.mesh,'faces','FormBasis','jumpgradwxn_S');
            val = FormBasis.sideValues('Grad',pts,pts_ref,parts,fe);
            val = FormBasis.stackSides(FormBasis.getJump(val));
            normals = fe.mesh.normal4f;
            normals = permute(normals(parts,:),[1 3 2]).*ones(1,size(val,2));
            val = cross(val,normals,3);
        end

        function val = jumpsymcurlwt_S(pts,pts_ref,parts,fe)
            %jumpsymcurlwt_S - Tangential jump of the symmetric curl
            %   VAL = jumpsymcurlwt_S(PTS,PTS_REF,PARTS,FE) returns the
            %   jump of the symmetric part of the curl of the 2d vector
            %   basis of FE (see symcurlw_T) times the tangents across the
            %   edges PARTS (edges only, see jump).
            val = FormBasis.jump('SymCurl','t',pts,pts_ref,parts,fe);
        end

        function val = jumpd2wt_S(pts,pts_ref,parts,fe)
            %jumpd2wt_S - Tangential jump of the Hessian
            %   VAL = jumpd2wt_S(PTS,PTS_REF,PARTS,FE) returns the jump of
            %   the Hessian of the scalar basis of FE times the tangents,
            %   [D^2 v]*t, across the edges PARTS (edges only, see jump).
            val = FormBasis.jump('D2','t',pts,pts_ref,parts,fe);
        end

        function val = jumpdivd2wt_S(pts,pts_ref,parts,fe)
            %jumpdivd2wt_S - Tangential jump of the gradient of the Laplacian
            %   VAL = jumpdivd2wt_S(PTS,PTS_REF,PARTS,FE) returns the jump
            %   of div(D^2 v)*t = grad(laplace(v))*t of the scalar basis of
            %   FE across the edges PARTS (edges only, see jump).
            val = FormBasis.jump('DivD2','t',pts,pts_ref,parts,fe);
        end

        function val = jumpdivw_S(pts,pts_ref,parts,fe)
            %jumpdivw_S - Jump of the divergence
            %   VAL = jumpdivw_S(PTS,PTS_REF,PARTS,FE) returns the
            %   difference of the divergence of the basis of FE on the two
            %   neighboring elements of the sides PARTS, stacked along
            %   dim 2 (see jump).
            val = FormBasis.jump('Div','',pts,pts_ref,parts,fe);
        end

        function val = jumpdivwn_S(pts,pts_ref,parts,fe)
            %jumpdivwn_S - Normal jump of the divergence
            %   VAL = jumpdivwn_S(PTS,PTS_REF,PARTS,FE) returns the jump of
            %   the divergence of the basis of FE times the normals across
            %   the sides PARTS (see jump).
            val = FormBasis.jump('Div','n',pts,pts_ref,parts,fe);
        end

        function val = jumpdivwt_S(pts,pts_ref,parts,fe)
            %jumpdivwt_S - Tangential jump of the divergence
            %   VAL = jumpdivwt_S(PTS,PTS_REF,PARTS,FE) returns the jump of
            %   the divergence of the basis of FE times the tangents across
            %   the edges PARTS (edges only, see jump).
            val = FormBasis.jump('Div','t',pts,pts_ref,parts,fe);
        end

        function val = averagew_S(pts,pts_ref,parts,fe)
            %averagew_S - Mean of the basis across the sides
            %   VAL = averagew_S(PTS,PTS_REF,PARTS,FE) returns the mean
            %   of the basis values of FE over the two neighboring
            %   elements of the sides PARTS, stacked along dim 2.
            val = FormBasis.average('','',pts,pts_ref,parts,fe);
        end

        function val = averagegradw_S(pts,pts_ref,parts,fe)
            %averagegradw_S - Mean of the basis gradient
            %   VAL = averagegradw_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of the gradient of the basis of FE over the two
            %   neighboring elements of the sides PARTS, stacked along
            %   dim 2.
            val = FormBasis.average('Grad','',pts,pts_ref,parts,fe);
        end

        function val = averagecurlw_S(pts,pts_ref,parts,fe)
            %averagecurlw_S - Mean of the curl across the sides
            %   VAL = averagecurlw_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of the curl of the basis of FE over the two
            %   neighboring elements of the sides PARTS, stacked along
            %   dim 2.
            val = FormBasis.average('Curl','',pts,pts_ref,parts,fe);
        end

        function val = averagewn_S(pts,pts_ref,parts,fe)
            %averagewn_S - Mean of the normal basis
            %   VAL = averagewn_S(PTS,PTS_REF,PARTS,FE) returns the mean
            %   of the basis of FE times the normals over the two
            %   neighboring elements of the sides PARTS (see average).
            val = FormBasis.average('','n',pts,pts_ref,parts,fe);
        end

        function val = averagewt_S(pts,pts_ref,parts,fe)
            %averagewt_S - Mean of the tangential basis
            %   VAL = averagewt_S(PTS,PTS_REF,PARTS,FE) returns the mean
            %   of the basis of FE times the tangents over the two
            %   neighboring elements of the edges PARTS (edges only,
            %   see average).
            val = FormBasis.average('','t',pts,pts_ref,parts,fe);
        end

        function val = averagegradwn_S(pts,pts_ref,parts,fe)
            %averagegradwn_S - Mean of the normal basis gradient
            %   VAL = averagegradwn_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of the normal derivative of the basis of FE over the
            %   two neighboring elements of the sides PARTS (see
            %   average).
            val = FormBasis.average('Grad','n',pts,pts_ref,parts,fe);
        end

        function val = averagegradwt_S(pts,pts_ref,parts,fe)
            %averagegradwt_S - Mean of the tangential basis gradient
            %   VAL = averagegradwt_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of the tangential derivative of the basis of FE over
            %   the two neighboring elements of the edges PARTS (edges
            %   only, see average).
            val = FormBasis.average('Grad','t',pts,pts_ref,parts,fe);
        end

        function val = averagewnn_S(pts,pts_ref,parts,fe)
            %averagewnn_S - Mean of the normal-normal component
            %   VAL = averagewnn_S(PTS,PTS_REF,PARTS,FE) returns the mean
            %   of the normal-normal components of the matrix-valued
            %   basis of FE over the two neighboring elements of the
            %   sides PARTS (see average).
            val = FormBasis.average('','nn',pts,pts_ref,parts,fe);
        end

        function val = averaged2wn_S(pts,pts_ref,parts,fe)
            %averaged2wn_S - Mean of the Hessian times the normal
            %   VAL = averaged2wn_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of the Hessian of the basis of FE times the
            %   normals, {D^2 w}*n, over the two neighboring elements of
            %   the sides PARTS (see average).
            val = FormBasis.average('D2','n',pts,pts_ref,parts,fe);
        end

        function val = averaged2wt_S(pts,pts_ref,parts,fe)
            %averaged2wt_S - Mean of the Hessian times the tangent
            %   VAL = averaged2wt_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of the Hessian of the basis of FE times the
            %   tangents, {D^2 w}*t, over the two neighboring elements of
            %   the edges PARTS (edges only, see average).
            val = FormBasis.average('D2','t',pts,pts_ref,parts,fe);
        end

        function val = averaged2wnn_S(pts,pts_ref,parts,fe)
            %averaged2wnn_S - Mean of the second normal derivative
            %   VAL = averaged2wnn_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of the second normal derivative of the basis of FE
            %   over the two neighboring elements of the sides PARTS
            %   (see average).
            val = FormBasis.average('D2','nn',pts,pts_ref,parts,fe);
        end

        function val = averagedivd2wn_S(pts,pts_ref,parts,fe)
            %averagedivd2wn_S - Mean of the normal derivative of the Laplacian
            %   VAL = averagedivd2wn_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of div(D^2 v)*n = grad(laplace(v))*n of the scalar
            %   basis of FE over the two neighboring elements of the
            %   sides PARTS (see average).
            val = FormBasis.average('DivD2','n',pts,pts_ref,parts,fe);
        end

        function val = averagedivd2wt_S(pts,pts_ref,parts,fe)
            %averagedivd2wt_S - Mean tangential derivative of the Laplacian
            %   VAL = averagedivd2wt_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of div(D^2 v)*t = grad(laplace(v))*t of the scalar
            %   basis of FE over the two neighboring elements of the
            %   edges PARTS (edges only, see average).
            val = FormBasis.average('DivD2','t',pts,pts_ref,parts,fe);
        end

        function val = averagedivwn_S(pts,pts_ref,parts,fe)
            %averagedivwn_S - Mean of the divergence times the normal
            %   VAL = averagedivwn_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of the divergence of the basis of FE times the
            %   normals over the two neighboring elements of the sides
            %   PARTS (see average).
            val = FormBasis.average('Div','n',pts,pts_ref,parts,fe);
        end

        function val = averagedivwt_S(pts,pts_ref,parts,fe)
            %averagedivwt_S - Mean of the divergence times the tangent
            %   VAL = averagedivwt_S(PTS,PTS_REF,PARTS,FE) returns the
            %   mean of the divergence of the basis of FE times the
            %   tangents over the two neighboring elements of the edges
            %   PARTS (edges only, see average).
            val = FormBasis.average('Div','t',pts,pts_ref,parts,fe);
        end

        %% side values
        function val = sideValues(name,pts,pts_ref,parts,fe)
            %sideValues - Basis values on both neighboring elements
            %   VAL = sideValues(NAME,PTS,PTS_REF,PARTS,FE) returns the
            %   basis values of FE on the two neighboring elements of the
            %   sides PARTS (see Form.sideOf), with the elements along
            %   the last dimension (zero for the missing second element
            %   of a boundary side). NAME selects the derivative: ''
            %   (values), 'Grad', 'D2', 'D3', 'Curl' or 'Div'. FE provides
            %   them directly by getEdge<NAME>Basis or getFace<NAME>Basis,
            %   or they are reconstructed from its element basis
            %   get<NAME>Basis. Two values are derived from these: 'DivD2',
            %   div(D^2 v) = grad(laplace(v)) of scalar bases from 'D3'
            %   (see divd2w_T), and 'SymCurl', the symmetric part of the
            %   row-wise curl of 2d vector bases from 'Grad' (see
            %   symcurlw_T).
            switch name
                case 'DivD2'
                    d3 = FormBasis.sideValues('D3',pts,pts_ref,parts,fe);
                    if size(d3,3) == 4
                        val = cat(3,sum(d3(:,:,[1 4],:),3),sum(d3(:,:,[2 3],:),3));
                    else
                        val = cat(3,sum(d3(:,:,[1 6 8],:),3),sum(d3(:,:,[4 2 9],:),3), ...
                            sum(d3(:,:,[5 7 3],:),3));
                    end
                    return
                case 'SymCurl'
                    grad = FormBasis.sideValues('Grad',pts,pts_ref,parts,fe);
                    if size(fe.mesh.c4n,2) ~= 2 || ndims(grad) ~= 4 || size(grad,3) ~= 4
                        error('fem4m:FormBasis:symCurlNotDefined', ...
                            'FormBasis: sym(curl(u)) is only defined for 2d vector fields');
                    end
                    % grad = [dx v1, dx v2, dy v1, dy v2], curl column by column
                    curl = cat(3,grad(:,:,3:4,:),-grad(:,:,1:2,:));
                    val = FormType.matrix2d.symPart(curl,3);
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
            sideMethod = ['get',side,name,'Basis'];
            elemMethod = ['get',name,'Basis'];
            if ismethod(fe,sideMethod)
                val = fe.(sideMethod)(pts,pts_ref,parts);
            elseif ~fe.isRefElemBased
                val = elemValues(pts,parts,@(x,ind)fe.(elemMethod)(x,[],ind),fe.mesh);
            else
                val = refValues(pts_ref,@(x)fe.(elemMethod)([],x,[]),fe.mesh);
                val = val(parts,:,:,:,:);
            end
        end

        %% side operators
        function val = oneSide(name,letters,pts,pts_ref,parts,fe)
            %oneSide - Side operator on the first neighboring element
            %   VAL = oneSide(NAME,LETTERS,PTS,PTS_REF,PARTS,FE)
            %   evaluates the basis values NAME of FE (see sideValues)
            %   on the first neighboring element of the sides PARTS and
            %   multiplies them by the directions LETTERS of that
            %   element: none for '', the normal 'n' or tangent 't'
            %   (see oneSideTimes), or for two letters AB the matrix
            %   values between them, B'*M*A (see contractOneSide).
            dirs = Form.directions(letters,pts_ref,parts,fe.mesh);
            dirs = cellfun(@Form.getT1,dirs,'UniformOutput',false);
            val = Form.getT1(FormBasis.sideValues(name,pts,pts_ref,parts,fe));
            if isscalar(dirs)
                val = FormBasis.oneSideTimes(val,dirs{1},fe);
            elseif numel(dirs) == 2
                val = FormBasis.contractOneSide(val,dirs{2},dirs{1},fe);
            end
        end

        function val = jump(name,letters,pts,pts_ref,parts,fe)
            %JUMP - Jump of a side operator across the sides
            %   VAL = JUMP(NAME,LETTERS,PTS,PTS_REF,PARTS,FE) returns the
            %   jump of the basis values NAME of FE (see sideValues)
            %   times the directions LETTERS (see oneSide) across the
            %   sides PARTS, with the basis functions of the second
            %   element following those of the first. Each element is
            %   multiplied by its own directions, which on the second
            %   element are the negated ones, or its conormals on
            %   surfaces. With one direction this already gives the
            %   jump; with none or two the signs cancel and the second
            %   element is negated.
            dirs = Form.directions(letters,pts_ref,parts,fe.mesh);
            val = FormBasis.sideValues(name,pts,pts_ref,parts,fe);
            if mod(numel(dirs),2) == 0
                val = FormBasis.getJump(val);
            end
            val = FormBasis.bothSidesTimes(val,dirs,fe);
        end

        function val = average(name,letters,pts,pts_ref,parts,fe)
            %AVERAGE - Mean of a side operator across the sides
            %   VAL = AVERAGE(NAME,LETTERS,PTS,PTS_REF,PARTS,FE) returns
            %   the mean of the basis values NAME of FE (see sideValues)
            %   over the two neighboring elements of the sides PARTS
            %   (see getAverage) times the directions LETTERS of the
            %   first element (see oneSide), with the basis functions of
            %   the second element following those of the first.
            e4s = Form.sideOf(fe.mesh,'e4s');
            dirs = Form.directions(letters,pts_ref,parts,fe.mesh);
            dirs = cellfun(@Form.getT1,dirs,'UniformOutput',false);
            val = FormBasis.sideValues(name,pts,pts_ref,parts,fe);
            val = FormBasis.getAverage(val,e4s(parts,:));
            val = FormBasis.bothSidesTimes(val,dirs,fe);
        end

        %% value type
        function type = typeOf(val,fe)
            %typeOf - Value type of element basis values
            %   TYPE = typeOf(VAL,FE) returns the FormType of the basis
            %   values VAL of FE on elements or on one side of the
            %   edges or faces, stored as rows x basis functions x
            %   components. Nonsingleton third and fourth dimensions
            %   give FormType.tensor, as does a number of components
            %   that matches no other type (see FormType.of).
            type = FormType.of(size(val,[3 4]),size(fe.mesh.c4n,2));
        end

        function type = sideTypeOf(val,fe)
            %sideTypeOf - Value type of two-sided basis values
            %   TYPE = sideTypeOf(VAL,FE) returns the FormType of the
            %   basis values VAL of FE on both sides of edges or faces,
            %   stored as rows x basis functions x components x 2
            %   sides; scalar values have no component dimension.
            sz = size(val);
            type = FormType.of(sz(3:end-1),size(fe.mesh.c4n,2));
        end

        %% one-sided contractions
        function val = oneSideTimes(basis,dirs,fe)
            %oneSideTimes - One-sided basis values times directions
            %   VAL = oneSideTimes(BASIS,DIRS,FE) multiplies the basis
            %   values BASIS of FE on one side (see typeOf) by the
            %   directions DIRS, rows x dim: scalars times DIRS, the
            %   inner product for vectors and the matrix times DIRS for
            %   matrices.
            type = FormBasis.typeOf(basis,fe);
            if type.isMatrix()
                val = sum(type.toMatrix(basis,3).*permute(dirs,[1 3 4 2]),4);
            elseif type.isVector()
                val = sum(basis.*permute(dirs,[1 3 2]),3);
            else
                val = basis.*permute(dirs,[1 3 2]);
            end
        end

        function val = contractOneSide(basis,left,right,fe)
            %contractOneSide - One-sided matrix values between two directions
            %   VAL = contractOneSide(BASIS,LEFT,RIGHT,FE) returns
            %   LEFT'*M*RIGHT for the matrix-valued basis values M in
            %   BASIS of FE on one side (see typeOf) and the directions
            %   LEFT and RIGHT, rows x dim.
            type = FormBasis.typeOf(basis,fe);
            weights = permute(left,[1 3 2]).*permute(right,[1 3 4 2]);
            val = sum(type.toMatrix(basis,3).*weights,[3 4]);
        end

        %% two-sided contractions
        function val = sideTimes(basis,dirs,fe)
            %sideTimes - Two-sided basis values times directions
            %   VAL = sideTimes(BASIS,DIRS,FE) multiplies the two-sided
            %   basis values BASIS of FE (see sideTypeOf) by the
            %   directions DIRS, rows x dim x 2 sides or rows x dim for
            %   both sides: scalars times DIRS, the inner product for
            %   vectors and the matrix times DIRS for matrices. The
            %   basis functions of the second side follow those of the
            %   first along dim 2.
            type = FormBasis.sideTypeOf(basis,fe);
            if type.isMatrix()
                val = sum(type.toMatrix(basis,3).*permute(dirs,[1 4 5 2 3]),4);
            elseif type.isVector()
                val = sum(basis.*permute(dirs,[1 4 2 3]),3);
            else
                val = permute(basis,[1 2 4 3]).*permute(dirs,[1 4 2 3]);
            end
            val = FormBasis.stackSides(val);
        end

        function val = contractSides(basis,left,right,fe)
            %contractSides - Two-sided matrix values between two directions
            %   VAL = contractSides(BASIS,LEFT,RIGHT,FE) returns
            %   LEFT'*M*RIGHT for the matrix-valued two-sided basis
            %   values M in BASIS of FE (see sideTypeOf) and the
            %   directions LEFT and RIGHT, rows x dim x 2 sides or
            %   rows x dim for both sides. The basis functions of the
            %   second side follow those of the first along dim 2.
            type = FormBasis.sideTypeOf(basis,fe);
            weights = permute(left,[1 4 2 5 3]).*permute(right,[1 4 5 2 3]);
            val = sum(type.toMatrix(basis,3).*weights,[3 4]);
            val = FormBasis.stackSides(val);
        end

        %% element rows
        function val = elementRows(val,parts,fe)
            %elementRows - Keep the rows of the elements PARTS
            %   VAL = elementRows(VAL,PARTS,FE) keeps the rows PARTS of the
            %   element values VAL of a reference element based FE, which
            %   may evaluate all elements and ignore PARTS. VAL is returned
            %   unchanged unless FE is reference element based and VAL has
            %   one row per element of the mesh.
            if fe.isRefElemBased && size(val,1) == fe.mesh.nrElems
                val = val(parts,:,:,:,:);
            end
        end

        %% jumps & averages (basis)
        function val = getJump(val)
            %getJump - Negate the second side values
            %   VAL = getJump(VAL) keeps both side values of VAL, stored
            %   along the last dimension, with the second one negated,
            %   so that stacking them (see stackSides) represents the
            %   jump. The size of VAL is kept.
            side2 = [repmat({':'},1,ndims(val)-1),{2}];
            val(side2{:}) = -val(side2{:});
        end

        function val = stackSides(val)
            %stackSides - Stack the two side values along the basis
            %   VAL = stackSides(VAL) turns the two-sided values VAL,
            %   rows x basis functions x components x 2 sides, into
            %   rows x 2*basis functions x components, with the basis
            %   functions of the second side following those of the
            %   first.
            sz = size(val);
            val = permute(val,[1 2 numel(sz) 3:numel(sz)-1]);
            val = reshape(val,[sz(1),2*sz(2),sz(3:end-1)]);
        end

        function val = getAverage(val,e4s)
            %getAverage - Average side values of the basis
            %   VAL = getAverage(VAL,E4S) halves both side values of
            %   VAL, stored along the last dimension, so that the
            %   stacked values (see stackSides) give the mean, and keeps
            %   the full values on boundary sides, where the second
            %   entry of E4S is zero. E4S holds the neighboring elements
            %   of the sides of VAL, one row per side. The size of VAL
            %   is kept.
            originDim = size(val);
            val = 0.5*val(:,:);
            % value on the boundary
            val(e4s(:,2)==0,:) = 2.*val(e4s(:,2)==0,:);
            val = reshape(val,originDim);
        end

    end

    methods (Static, Access = private)

        function val = bothSidesTimes(val,dirs,fe)
            %bothSidesTimes - Two-sided values times their directions
            %   VAL = bothSidesTimes(VAL,DIRS,FE) stacks the two-sided
            %   basis values VAL of FE without directions (see
            %   stackSides), multiplies them by one direction (see
            %   sideTimes) or returns the matrix values between two
            %   directions, DIRS{2}'*M*DIRS{1} (see contractSides).
            if isempty(dirs)
                val = FormBasis.stackSides(val);
            elseif isscalar(dirs)
                val = FormBasis.sideTimes(val,dirs{1},fe);
            else
                val = FormBasis.contractSides(val,dirs{2},dirs{1},fe);
            end
        end

    end
end
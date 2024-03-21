classdef AWk < AW1
    %AWk - Arnold-Winther element of degree k for elasticity
    %   OBJ = AWk(DEGREE,MESH) creates the Arnold-Winther symmetric
    %   H(div)-conforming stress finite element of polynomial degree
    %   DEGREE >= 3 on the triangular mesh MESH, built from BDM functions
    %   whose divergence degree is at most k = DEGREE-2, the degree of
    %   the displacements. It has (3*k^2+17*k+28)/2 degrees of freedom
    %   per element: the values of sigma at the vertices (3 per
    %   vertex), the moments of degree 0 to k of sigma*n on every edge
    %   (2+2*k per edge), and (3*k^2+5*k-2)/2 element moments (see
    %   prepareCoefficients). AWk(3,MESH) is the space of AW1.
    %
    %   AWk functions:
    %       prepareCoefficients - Compute the basis expansion coefficients
    %       getBasis            - Evaluate the Arnold-Winther basis functions
    %       getDivBasis         - Evaluate the divergence of the basis
    %       moment              - Evaluate the edge moment integrands
    %       momentG             - Evaluate the Neumann moment integrand
    %       boundaryMoments     - Assemble the Neumann boundary moments
    %
    %   AWk properties:
    %       bdm           - BDM matrix element used for the construction
    %       k             - Polynomial degree of the moments, DEGREE-2
    %       divkIndx      - Logical index of BDM dofs with div degree <= k
    %       divkBasisFunc - Number of selected BDM basis functions
    %
    %   See also AW1, AWH, AW0, HZ1, BDMk, DGDk

    % Copyright (c) 2026 by the fem4m authors

    properties
        bdm           % BDM matrix element used for the construction
        k             % polynomial degree of the moments, DEGREE-2
        divkIndx      % logical index of BDM dofs with div degree <= k
        divkBasisFunc % number of selected BDM basis functions
    end

    methods

        function obj = AWk(degree,mesh)
            %AWk - Create an Arnold-Winther finite element object
            if nargin < 2 || ~isscalar(degree) || degree < 3
                error(['fem4m:',class(obj),':invalidDegree'], ...
                    '%s degree must be >= 3, the stresses of degree DEGREE have divergence of degree DEGREE-2', ...
                    class(obj));
            end
            obj.degree = degree;
            obj.k = degree-2;
            obj.mesh = mesh;
            obj.nrBasisFunc = 9 + 6 + 6*obj.k + (3*obj.k^2 + 5*obj.k - 2)/2;
            obj.bdm = BDMk(obj.degree,obj.mesh);
            indx = obj.bdm.divDegree4DoF <= obj.k;
            obj.divkIndx = [indx,indx];
            obj.divkBasisFunc = nnz(obj.divkIndx);
            obj.bdm = FEMatrix(obj.bdm,obj.bdm);
            obj = obj.prepareCoefficients();
            obj = obj.dofCounts();
            obj = obj.dofManager();
            obj = obj.boundaryConstraints();
        end

        %% Arnold-Winther transformation from BDM functions with div <= k
        function obj = prepareCoefficients(obj)
            %prepareCoefficients - Compute the basis expansion coefficients
            %   OBJ = prepareCoefficients(OBJ) solves on every element the
            %   local Vandermonde system given by the nodal values, the
            %   edge moments of sigma*n of degree 0 to k, the element
            %   moments (the means of sigma and, for k >= 2, the moments
            %   against eps(v) for the vector fields v whose components
            %   have degree 2 to k, and against curlcurl of the B2k
            %   bubbles), and the symmetry constraints (the skew part of
            %   sigma vanishes) to obtain the expansion of the AWk basis
            %   in the selected BDM basis (basisCoefficients).
            % 9 nodal dofs
            nodes = zeros(obj.mesh.nrElems,obj.divkBasisFunc,9);
            z = [0,0;1,0;0,1];
            for j = 1 : 3
                val = obj.bdm.getBasis([],z(j,:),[]);
                nodes(:,:,3*(j-1)+(1:2)) = val(:,obj.divkIndx,[1 4]);
                nodes(:,:,3*(j-1)+3) = 1/2*(val(:,obj.divkIndx,2)+val(:,obj.divkIndx,3));
            end
            % edge moments 0,1,...,k
            edgeMoments = integrate(@(~,refp,~)obj.moment(refp,obj.k),...
                obj.mesh.n4ed,1:obj.mesh.nrEdges,obj.mesh,obj.degree+obj.k)./obj.mesh.length4ed;
            % elem moments 0
            dg = DGDk(0,obj.mesh);
            dg = FEMatrix(dg,dg,dg,dg);
            form = BilinearForm('u:vdx');
            [~,elemMoments] = form.assemble(dg,obj.bdm);
            % symmetry constraints
            skewdg = FESkewSymMatrix(DGDk(obj.degree,obj.mesh));
            [~,symConstraints] = form.assemble(skewdg,obj.bdm);
            nodes = permute(nodes,[3 2 1]);
            elemMoments = permute(elemMoments./obj.mesh.area4e,[2 3 1]);
            symConstraints = permute(symConstraints./obj.mesh.area4e,[2 3 1]);
            % higher order elem moments
            if obj.k >= 2
                dgeps = DGDk(obj.k,obj.mesh);
                dgeps = FEVector(dgeps,dgeps);
                form = BilinearForm('eps(u):v dx');
                [~,elemMomentsEps] = form.assemble(dgeps,obj.bdm);
                elemMomentsEps = permute(elemMomentsEps./sqrt(obj.mesh.area4e),[2 3 1]);
                dgCurlCurl = B2k(obj.k+4,obj.mesh);
                form = BilinearForm('curlcurl(u):v dx');
                [~,elemMomentsCurlCurl] = form.assemble(dgCurlCurl,obj.bdm);
                elemMomentsCurlCurl = permute(elemMomentsCurlCurl,[2 3 1]);
            end
            obj.basisCoefficients = zeros(obj.mesh.nrElems,obj.nrBasisFunc,obj.divkBasisFunc);
            for curElem = 1 : obj.mesh.nrElems
                C = zeros(obj.divkBasisFunc,obj.divkBasisFunc);
                % 9 node
                C(1:9,:) = nodes(1:9,:,curElem);
                % 6+6*obj.k edges
                TplusMinus = (obj.mesh.e4ed(obj.mesh.ed4e(curElem,:),:)==curElem);
                for j = 0 : obj.k
                    for m = 1 : 3
                        C(9+1+j*2+2*(obj.k+1)*(m-1),:) = edgeMoments(obj.mesh.ed4e(curElem,m),:,1,TplusMinus(m,:),j+1);
                        C(9+1+j*2+2*(obj.k+1)*(m-1)+1,:) = edgeMoments(obj.mesh.ed4e(curElem,m),:,2,TplusMinus(m,:),j+1);
                    end
                end
                % 3 moments order 0 -> eps(P_1^2) without RBM
                C(9+6+6*obj.k+(1:3),:) = 1/2*(elemMoments([1 4 2],obj.divkIndx,curElem)+elemMoments([1 4 3],obj.divkIndx,curElem));
                if obj.k>=2
                    epsIndx = dgeps.fe{1}.degree4DoF>1;
                    C(9+6+6*obj.k+3+(1:dgeps.nrBasisFunc-6),:) = elemMomentsEps([epsIndx,epsIndx],obj.divkIndx,curElem);
                    C(9+6+6*obj.k+3+dgeps.nrBasisFunc-6+(1:dgCurlCurl.nrBasisFunc),:) = elemMomentsCurlCurl(:,obj.divkIndx,curElem);
                end
                % symmetry constraints
                C(end-skewdg.nrBasisFunc+1:end,:) = symConstraints(:,obj.divkIndx,curElem);
                % compute coefficients
                x = C\[eye(obj.nrBasisFunc);zeros(obj.divkBasisFunc-obj.nrBasisFunc,obj.nrBasisFunc)];
                obj.basisCoefficients(curElem,:,:) = x';
            end

        end

        %% Stress basis functions
        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the Arnold-Winther basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the symmetric
            %   Arnold-Winther stress basis functions at the point X
            %   (1-by-2) of the reference triangle in all elements. VAL
            %   is nrElems-by-nrBasisFunc-by-3 with the components
            %   [s11 s22 s12] (Voigt notation), each basis function
            %   divided by the scaling of its degree of freedom. PTS and
            %   PARTS are ignored.
            basis = obj.bdm.getBasis([],x,[]);
            basis = basis(:,obj.divkIndx,:);
            c = obj.basisCoefficients;
            val = matMul(c,basis);
            val = 1/2*(val(:,:,[1 4 2])+val(:,:,[1 4 3]));
            val = val./obj.scaling(obj.dof4e);
        end

        function val = getDivBasis(obj,~,x,~)
            %getDivBasis - Evaluate the divergence of the basis
            %   VAL = getDivBasis(OBJ,PTS,X,PARTS) evaluates the row-wise
            %   divergences of the Arnold-Winther stress basis functions
            %   at the point X (1-by-2) of the reference triangle in all
            %   elements. VAL is nrElems-by-nrBasisFunc-by-2, scaled as
            %   in getBasis. PTS and PARTS are ignored.
            basis = obj.bdm.getDivBasis([],x,[]);
            basis = basis(:,obj.divkIndx,:);
            c = obj.basisCoefficients;
            val = matMul(c,basis);
            val = val./obj.scaling(obj.dof4e);
        end

        %% Moments
        function val = moment(obj,refp,p)
            %MOMENT - Evaluate the edge moment integrands
            %   VAL = MOMENT(OBJ,REFP,P) evaluates, at the edge reference
            %   point REFP, the traction sigma*n (with the edge normals
            %   normal4ed) of the symmetric part sigma of the selected BDM
            %   basis functions times the Legendre polynomials of degree
            %   0 to P, the integrands of the edge moments. VAL is
            %   nrEdges-by-divkBasisFunc-by-2-by-2-by-(P+1): the two
            %   traction components, the two elements e4ed of each edge,
            %   and the degree of the moment.
            sigma = getEdgeValues4RefElemFunktion(refp,@(z)obj.bdm.getBasis([],z,[]),obj.mesh);
            sigma = sigma(:,obj.divkIndx,:,:);
            n = obj.mesh.normal4ed;
            val(:,:,1,:) = sigma(:,:,1,:).*n(:,1) + 1/2*(sigma(:,:,2,:)+sigma(:,:,3,:)).*n(:,2);
            val(:,:,2,:) = 1/2*(sigma(:,:,3,:)+sigma(:,:,2,:)).*n(:,1) + sigma(:,:,4,:).*n(:,2);
            val2 = reshape(JacobiPolynomials.poly((2*refp(1)-1),p,0,0,0),[1,1,1,1,p+1]);
            val = val.*val2;
        end

        function val = momentG(obj,pts,refp,ind,g,p)
            %momentG - Evaluate the Neumann moment integrand
            %   VAL = momentG(OBJ,PTS,REFP,IND,G,P) evaluates the
            %   traction G(PTS,NORMALS) on the edges IND times the
            %   Legendre polynomial of degree P at the edge reference
            %   point REFP, the integrand of the edge moment of degree P
            %   of the Neumann data.
            normals = obj.mesh.normal4ed;
            normals = normals(ind,:);
            curG = g(pts,normals);
            val2 = JacobiPolynomials.poly((2*refp(1)-1),p,0,0,0);
            val = curG.*val2(end);
        end

        function val = boundaryMoments(obj,g,gDegree)
            %boundaryMoments - Assemble the Neumann boundary moments
            %   VAL = boundaryMoments(OBJ,G) computes the edge moments of
            %   degree 0 to k of the traction function handle
            %   G(PTS,NORMALS) on all Neumann boundary edges and returns
            %   them at the Neumann edge dofs of the vector VAL
            %   (nrDoF-by-1, zero elsewhere).
            %
            %   VAL = boundaryMoments(OBJ,G,GDEGREE) also specifies the
            %   polynomial degree GDEGREE of the data G, which defaults
            %   to OBJ.degree and only sets the quadrature degree: the
            %   moment of order j-1 is integrated with the degree
            %   max(momentDegree,j+GDEGREE), at least with the Gauss rule
            %   of AW1, so a GDEGREE below momentDegree-j has no
            %   effect. AWk(3,MESH) and AWH(3,MESH) thus give the same
            %   boundary function as AW1; their first order moment dofs
            %   are twice those of AW1.
            if nargin < 3 || isempty(gDegree)
                gDegree = obj.degree;
            end
            val = zeros(obj.nrDoF,1);
            ed4Nb = obj.mesh.ed4Nb;
            length4ed = obj.mesh.length4ed;
            if ~isempty(obj.mesh.n4Nb)
                for j = 1 : obj.k+1
                    moments = (1./length4ed(obj.mesh.ed4Nb)*ones(1,2)).*integrate(...
                        @(points,refp,ind)obj.momentG(points,refp,ind,g,j-1),obj.mesh.n4Nb,obj.mesh.ed4Nb,obj.mesh, ...
                        max(obj.momentDegree,j+gDegree));
                    fixNodes0 = 3*obj.mesh.nrNodes+2*obj.mesh.nrEdges*(j-1)+[ed4Nb,obj.mesh.nrEdges+ed4Nb];
                    val(fixNodes0(:)) = moments(:);
                end
            end

        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of the 3 nodal stress values,
            %   2+2*k edge moments, and the element moments.
            obj.nrBasis4n = 3;
            obj.nrBasis4ed = 2+2*obj.k;
            obj.nrBasis4e = obj.nrBasisFunc-9-3*(2+2*obj.k);
        end

    end
end
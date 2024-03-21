classdef ProblemTest < UnitTestCase
    %ProblemTest - Consistency tests of the problem definitions
    %   The derivatives of the exact solutions of the problems agree with
    %   central differences of the functions they derive from, and the
    %   right-hand sides are the differential operators of the problems
    %   applied to the exact solutions, evaluated with central differences
    %   of the exact derivatives at interior points of the domains:
    %
    %       scalar elliptic   f = -div(alpha grad u) + beta.grad u + gamma u
    %                         (outside the contact set of obstacle problems)
    %       vector elliptic   f = -Laplace u
    %       biharmonic        f = Laplace^2 u
    %       elasticity        f = -div sigma, sigma = C eps(u)
    %       Stokes            f = -nu Laplace u + grad p, div u = 0
    %       Helmholtz         f = -Laplace u - omega^2 n^2 u
    %       Maxwell (2d)      f = curl(mu^-1 rot u) + alpha u
    %       Maxwell (3d)      f = curl(mu^-1 curl u) + alpha u
    %       surface           f = -Laplace_Gamma u on the surface
    %
    %   In 2d, curl w = (d_2 w, -d_1 w) for scalar fields w and
    %   rot v = d_1 v_2 - d_2 v_1 for vector fields v, row-wise for
    %   matrix fields. The row-wise rots of the strain and stress of the
    %   elasticity and Stokes problems and the rot of the source of the
    %   Maxwell problems are checked as well.
    %
    %   The paraxial Gaussian beam of problemHelmholtzGaussianBeam is only
    %   an approximate solution, with a gradient that agrees with central
    %   differences and a relative Helmholtz residual below one percent.
    %
    %   The surfaces are zero sets of their surface functions with
    %   consistent gradients and Hessians. The reference eigenvalues are
    %   positive and nondecreasing, the eigenvalues of several modes are
    %   those of the single modes, and the problems without an exact
    %   solution give finite data of the right size.
    %
    %   The twisted diffusion eigenvalue problem has a symmetric diffusion
    %   matrix alpha = inv(DPhi)*inv(DPhi)' with det(alpha) = det(DPhi) = 1,
    %   the derivatives dalpha, DPhi and exactGradU agree with central
    %   differences, the twist fixes the boundary, and the eigenfunctions
    %   of the eigenspaces satisfy -div(alpha grad u) = lambda u.
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   ProblemBiharmonicTemplateExact, ProblemElasticityTemplateExact,
    %   ProblemStokesTemplateExact, ProblemHelmholtzTemplateExact,
    %   ProblemSurfaceTemplate, ProblemEigenvalue, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        scalarElliptic = struct( ...
            'CircleExact',@()problemScalarEllipticCircleExact, ...
            'CubeExact',@()problemScalarEllipticCubeExact, ...
            'CubeExact2',@()problemScalarEllipticCubeExact2, ...
            'CubeMonomExact',@()problemScalarEllipticCubeMonomExact(2,1,1), ...
            'LbrickExact',@()problemScalarEllipticLbrickExact, ...
            'LshapeExact',@()problemScalarEllipticLshapeExact, ...
            'LshapeOscExact',@()problemScalarEllipticLshapeOscExact, ...
            'SlitExact',@()problemScalarEllipticSlitExact, ...
            'SquareClampedExact',@()problemScalarEllipticSquareClampedExact, ...
            'SquareExact',@()problemScalarEllipticSquareExact, ...
            'SquareFullEllipticExact',@()problemScalarEllipticSquareFullEllipticExact, ...
            'SquareKelloggExact',@()problemScalarEllipticSquareKelloggExact, ...
            'SquareMonomExact',@()problemScalarEllipticSquareMonomExact(2,1), ...
            'SquareMountainValleyExact',@()problemScalarEllipticSquareMountainValleyExact, ...
            'SquareNeumannExact',@()problemScalarEllipticSquareNeumannExact, ...
            'SquareOscExact',@()problemScalarEllipticSquareOscExact, ...
            'SquareWaterfallExact',@()problemScalarEllipticSquareWaterfallExact, ...
            'ObstacleLshapeExact',@()problemObstacleLshapeExact, ...
            'ObstacleSquareExact',@()problemObstacleSquareExact, ...
            'ObstacleSquareExactSmooth',@()problemObstacleSquareExactSmooth, ...
            'ConvectionDominatedSquareExact',@()problemConvectionDominatedSquareExact, ...
            'TransportSquareExact',@()problemTransportSquareExact, ...
            'TransportSquareMonomExact',@()problemTransportSquareMonomExact(2,1), ...
            'VectorEllipticSquareExact',@()problemVectorEllipticSquareExact)
        biharmonic = struct( ...
            'CubeClampedExact',@()problemBiharmonicCubeClampedExact, ...
            'CubeClampedExact2',@()problemBiharmonicCubeClampedExact2, ...
            'LshapeClampedExact',@()problemBiharmonicLshapeClampedExact, ...
            'SimplySupportedSingularRHSCubeExact',@()problemBiharmonicSimplySupportedSingularRHSCubeExact, ...
            'SquareClampedExact',@()problemBiharmonicSquareClampedExact, ...
            'SquareClampedExact2',@()problemBiharmonicSquareClampedExact2, ...
            'SquareClampedOscExact',@()problemBiharmonicSquareClampedOscExact, ...
            'SquareSimplySupportedExact',@()problemBiharmonicSquareSimplySupportedExact)
        elasticity = struct( ...
            'CircularInclusion',@()problemElasticityCircularInclusion('E',1,'nu',0.3), ...
            'LshapeExact',@()problemElasticityLshapeExact('E',1,'nu',0.3), ...
            'SquareExact',@()problemElasticitySquareExact('E',1,'nu',0.3), ...
            'SquareExactLocking',@()problemElasticitySquareExactLocking('E',1,'nu',0.3), ...
            'SquareExactZeroRhs',@()problemElasticitySquareExactZeroRhs('E',1,'nu',0.3), ...
            'SquareMonomExact',@()problemElasticitySquareMonomExact('E',1,'nu',0.3,2,1), ...
            'SquareNeumannExact',@()problemElasticitySquareNeumannExact('E',1,'nu',0.3))
        stokes = struct( ...
            'ExactCollidingFlow',@()problemStokesExactCollidingFlow, ...
            'ExactSlit',@()problemStokesExactSlit, ...
            'ExactSquareP3Sigma',@()problemStokesExactSquareP3Sigma, ...
            'ExactSquareSinus',@()problemStokesExactSquareSinus, ...
            'SquareMonomExact',@()problemStokesSquareMonomExact(2,2,1,1))
        helmholtz = struct( ...
            'LshapeExact',@()problemHelmholtzLshapeExact(2), ...
            'LshapeSmoothExact',@()problemHelmholtzLshapeSmoothExact(2), ...
            'PlaneWave',@()problemHelmholtzPlaneWave(2), ...
            'SquareExact',@()problemHelmholtzSquareExact(2), ...
            'SquareExact2',@()problemHelmholtzSquareExact2(2), ...
            'SquareMonomExact',@()problemHelmholtzSquareMonomExact(2,2,1), ...
            'SquareReflection',@()problemHelmholtzSquareReflection(2))
        maxwell = struct( ...
            'CubeExact',@()problemMaxwellCubeExact, ...
            'CubeExact2',@()problemMaxwellCubeExact2, ...
            'ExactDblConnected',@()problemMaxwellExactDblConnected(2), ...
            'ExactLshape',@()problemMaxwellExactLshape(2), ...
            'ExactSquare',@()problemMaxwellExactSquare(2), ...
            'Lshape3dExact',@()problemMaxwellLshape3dExact, ...
            'SquareMonomExact',@()problemMaxwellSquareMonomExact(2,2,1))
        beamFrequency = struct('omega2',2,'omega8',8,'omega32',32)
        surfaceExact = struct( ...
            'EllipsoidExact',@()problemSurfaceEllipsoidExact, ...
            'SphereExact',@()problemSurfaceSphereExact, ...
            'SphereExact2',@()problemSurfaceSphereExact2, ...
            'SphereExact3',@()problemSurfaceSphereExact3, ...
            'SphereLExact',@()problemSurfaceSphereLExact, ...
            'SquareExact',@()problemSurfaceSquareExact, ...
            'TorusExact',@()problemSurfaceTorusExact)
        surface = struct( ...
            'Ellipsoid',@()problemSurfaceEllipsoid, ...
            'G88',@()problemSurfaceG88, ...
            'Sphere',@()problemSurfaceSphere, ...
            'SphereL',@()problemSurfaceSphereL, ...
            'Square',@()problemSurfaceSquare, ...
            'Torus',@()problemSurfaceTorus)
        eigenvalue = struct( ...
            'BucklingLshapeClamped',@()problemBiharmonicEigenvalueBucklingLshapeClamped, ...
            'BucklingSlitClamped',@()problemBiharmonicEigenvalueBucklingSlitClamped, ...
            'BucklingSquareClamped',@()problemBiharmonicEigenvalueBucklingSquareClamped, ...
            'VibrationLshapeClamped',@()problemBiharmonicEigenvalueVibrationLshapeClamped, ...
            'VibrationSlitClamped',@()problemBiharmonicEigenvalueVibrationSlitClamped, ...
            'VibrationSquareClamped',@()problemBiharmonicEigenvalueVibrationSquareClamped, ...
            'ElasticityCircle',@()problemElasticityEigenvalueCircle('lambda',1,'mu',1), ...
            'ElasticityCircleNeumann',@()problemElasticityEigenvalueCircleNeumann('lambda',1,'mu',1), ...
            'ElasticityCircleSb',@()problemElasticityEigenvalueCircleSb('lambda',1,'mu',1), ...
            'ElasticityCooks',@()problemElasticityEigenvalueCooks('lambda',1,'mu',1), ...
            'ElasticityLshape',@()problemElasticityEigenvalueLshape('lambda',1,'mu',1), ...
            'ElasticityLshapeNeumann',@()problemElasticityEigenvalueLshapeNeumann('lambda',1,'mu',1), ...
            'ElasticityLshapeSb',@()problemElasticityEigenvalueLshapeSb('lambda',1,'mu',1), ...
            'ElasticitySquare',@()problemElasticityEigenvalueSquare('lambda',1,'mu',1), ...
            'ElasticitySquareNeumann',@()problemElasticityEigenvalueSquareNeumann('lambda',1,'mu',1), ...
            'ElasticitySquareSb',@()problemElasticityEigenvalueSquareSb('lambda',1,'mu',1), ...
            'ElasticityStripeX',@()problemElasticityEigenvalueStripeX('lambda',1,'mu',1), ...
            'ElasticityStripeY',@()problemElasticityEigenvalueStripeY('lambda',1,'mu',1), ...
            'ScalarDblConnected',@()problemScalarEllipticEigenvalueDblConnected, ...
            'ScalarDumbbell',@()problemScalarEllipticEigenvalueDumbbell, ...
            'ScalarIsospectralDrumsA',@()problemScalarEllipticEigenvalueIsospectralDrumsA, ...
            'ScalarIsospectralDrumsB',@()problemScalarEllipticEigenvalueIsospectralDrumsB, ...
            'ScalarLshape',@()problemScalarEllipticEigenvalueLshape, ...
            'ScalarLshapeNeumann',@()problemScalarEllipticEigenvalueLshapeNeumann, ...
            'ScalarSlit',@()problemScalarEllipticEigenvalueSlit, ...
            'ScalarSquare',@()problemScalarEllipticEigenvalueSquare, ...
            'ScalarSquareNeumann',@()problemScalarEllipticEigenvalueSquareNeumann, ...
            'ScalarSquareTwist',@()problemScalarEllipticEigenvalueSquareTwist, ...
            'SurfaceEllipsoid',@()problemSurfaceEigenvalueEllipsoid, ...
            'SurfaceG88',@()problemSurfaceEigenvalueG88, ...
            'SurfaceSphere',@()problemSurfaceEigenvalueSphere, ...
            'SurfaceSphereL',@()problemSurfaceEigenvalueSphereL, ...
            'SurfaceSquare',@()problemSurfaceEigenvalueSquare, ...
            'SurfaceTorus',@()problemSurfaceEigenvalueTorus)
        data = struct( ...
            'BendingSquare',@()problemBendingSquare, ...
            'BiharmonicLshapeClamped',@()problemBiharmonicLshapeClamped, ...
            'ElasticityCooks',@()problemElasticityCooks('E',1,'nu',0.3), ...
            'ElasticityCrack',@()problemElasticityCrack('E',1,'nu',0.3), ...
            'MaxwellCube',@()problemMaxwellCube, ...
            'ObstacleSquareCusp',@()problemObstacleSquareCusp, ...
            'ObstacleSquareDist',@()problemObstacleSquareDist, ...
            'ScalarEllipticCube',@()problemScalarEllipticCube, ...
            'ScalarEllipticLbrick',@()problemScalarEllipticLbrick, ...
            'ScalarEllipticLshape',@()problemScalarEllipticLshape, ...
            'ScalarEllipticSquare',@()problemScalarEllipticSquare, ...
            'ScalarEllipticSquareMixedBd',@()problemScalarEllipticSquareMixedBd, ...
            'StokesBackwardFacingStep',@()problemStokesBackwardFacingStep, ...
            'StokesDrivenCavity',@()problemStokesDrivenCavity, ...
            'SurfaceStokesG88',@()problemSurfaceStokesG88, ...
            'SurfaceStokesSphereExact',@()problemSurfaceStokesSphereExact, ...
            'SurfaceStokesSphereExact2',@()problemSurfaceStokesSphereExact2, ...
            'SurfaceStokesSquareExact',@()problemSurfaceStokesSquareExact, ...
            'SurfaceVectorLaplaceSphereExact',@()problemSurfaceVectorLaplaceSphereExact, ...
            'SurfaceVectorLaplaceSphereExact2',@()problemSurfaceVectorLaplaceSphereExact2, ...
            'TransportSquare',@()problemTransportSquare, ...
            'TransportSquare2',@()problemTransportSquare2)
    end

    properties (Constant)
        h = 1e-5       % step of the central differences of first derivatives
        h2 = 1e-3      % step of the second differences of the biharmonic
        fdTol = 1e-6   % relative tolerance of the central differences
    end

    methods (Test)

        function testScalarElliptic(testCase,scalarElliptic)
            p = scalarElliptic();
            % the obstacle problems need fine points outside the contact set
            x = ProblemTest.interiorPoints(p.geom,1 + 2*isa(p,'ProblemObstacleTemplate'));
            u = p.exactU(x);

            ProblemTest.verifyFD(testCase,p.exactGradU(x),@p.exactU,x,'exactGradU');
            if size(u,2) == 1
                div = ProblemTest.divergence(@p.exactSigma,x);
                f = -div + sum(p.beta(x).*p.exactGradU(x),2) + p.gamma(x).*u;
            else
                f = -ProblemTest.laplace(@p.exactGradU,x);
            end
            if isa(p,'ProblemObstacleTemplate')
                % the equation holds outside the contact set, and the
                % multiplier L(u) - f is nonnegative at the lower and
                % nonpositive at the upper obstacle
                lowerContact = u <= p.lower(x) + 1e-6;
                upperContact = u >= p.upper(x) - 1e-6;
                multiplier = f - p.f(x);
                scale = max(1,max(abs(f)));
                testCase.verifyTrue(all(multiplier(lowerContact) >= -testCase.fdTol*scale), ...
                    'multiplier at the lower obstacle');
                testCase.verifyTrue(all(multiplier(upperContact) <= testCase.fdTol*scale), ...
                    'multiplier at the upper obstacle');
                free = ~lowerContact & ~upperContact;
                x = x(free,:);
                f = f(free,:);
            end
            if ~isempty(x)
                ProblemTest.verifyRelative(testCase,p.f(x),f,'f');
            end
        end

        function testBiharmonic(testCase,biharmonic)
            p = biharmonic();
            x = ProblemTest.interiorPoints(p.geom);
            d = size(x,2);

            ProblemTest.verifyFD(testCase,p.exactGradU(x),@p.exactU,x,'exactGradU');
            J = ProblemTest.jacobian(@p.exactGradU,x,testCase.h);
            if d == 2
                expected = [J(:,1,1),J(:,2,2),J(:,1,2)];
            else
                expected = [J(:,1,1),J(:,2,2),J(:,3,3),J(:,1,2),J(:,1,3),J(:,2,3)];
            end
            ProblemTest.verifyRelative(testCase,p.exactD2(x),expected,'exactD2');
            if ~contains(class(p),'SingularRHS')
                ProblemTest.verifyRelative(testCase,p.f(x), ...
                    ProblemTest.secondDifferences(@p.exactLaplaceU,x,testCase.h2),'f',1e-4);
            end
        end

        function testElasticity(testCase,elasticity)
            p = elasticity();
            x = ProblemTest.interiorPoints(p.geom);

            gradU = p.exactGradU(x);
            ProblemTest.verifyFD(testCase,gradU,@p.exactU,x,'exactGradU');
            eps = [gradU(:,1,1),gradU(:,2,2),(gradU(:,1,2) + gradU(:,2,1))/2];
            ProblemTest.verifyRelative(testCase,p.exactSigma(x),eps*p.C','exactSigma = C eps(u)');
            J = ProblemTest.jacobian(@p.exactSigma,x,testCase.h);  % [11 22 12] x derivative
            divSigma = [J(:,1,1) + J(:,3,2),J(:,3,1) + J(:,2,2)];
            ProblemTest.verifyRelative(testCase,p.f(x),-divSigma,'f');
            if ProblemTest.defines(p,'exactD2U')
                J = ProblemTest.jacobian(@(z)reshape(p.exactGradU(z),size(z,1),[]),x,testCase.h);
                J = reshape(J,[],2,2,2);  % component, derivative, derivative
                expected = cat(3,J(:,:,1,1),J(:,:,2,2),J(:,:,1,2));
                ProblemTest.verifyRelative(testCase,p.exactD2U(x),expected,'exactD2U');
            end
            ProblemTest.verifyRotEpsilon(testCase,p,x);
            if ismethod(p,'exactRotSigma')
                ProblemTest.verifyRelative(testCase,p.exactRotSigma(x), ...
                    ProblemTest.rot(@p.exactSigma,x),'exactRotSigma');
            end
            if ismethod(p,'exactRotRotSigma')
                ProblemTest.verifyRelative(testCase,p.exactRotRotSigma(x), ...
                    ProblemTest.rot(@p.exactRotSigma,x),'exactRotRotSigma');
            end
        end

        function testStokes(testCase,stokes)
            p = stokes();
            x = ProblemTest.interiorPoints(p.geom);

            gradU = p.exactGradU(x);
            ProblemTest.verifyFD(testCase,gradU,@p.exactU,x,'exactGradU');
            ProblemTest.verifyRelative(testCase,p.exactDivU(x),zeros(size(x,1),1),'div u');
            gradP = reshape(ProblemTest.jacobian(@p.exactP,x,testCase.h),size(x,1),[]);
            f = -p.nu*ProblemTest.laplace(@p.exactGradU,x) + gradP;
            ProblemTest.verifyRelative(testCase,p.f(x),f,'f');
            ProblemTest.verifyRotEpsilon(testCase,p,x);
        end

        function testHelmholtz(testCase,helmholtz)
            p = helmholtz();
            x = ProblemTest.interiorPoints(p.geom);

            ProblemTest.verifyFD(testCase,p.exactGradU(x),@p.exactU,x,'exactGradU');
            f = -ProblemTest.laplace(@p.exactGradU,x) ...
                - p.omega^2*p.refraction(x).^2.*p.exactU(x);
            ProblemTest.verifyRelative(testCase,p.f(x),f,'f');
        end

        function testGaussianBeam(testCase,beamFrequency)
            p = problemHelmholtzGaussianBeam(beamFrequency);
            x = ProblemTest.interiorPoints(p.geom);
            u = p.exactU(x);

            ProblemTest.verifyFD(testCase,p.exactGradU(x),@p.exactU,x,'exactGradU');
            residual = -ProblemTest.laplace(@p.exactGradU,x) - p.omega^2*u;
            testCase.verifyLessThan(max(abs(residual))/max(abs(p.omega^2*u)),1e-2, ...
                'paraxial approximation');
        end

        function testMaxwell(testCase,maxwell)
            p = maxwell();
            x = ProblemTest.interiorPoints(p.geom);
            if ismethod(p,'mu')
                mu = @(z)p.mu(z);
            else
                mu = @(z)ones(size(z,1),1);
            end
            if isprop(p,'alpha')
                alpha = p.alpha;
            else
                alpha = 0;
            end

            if size(x,2) == 3
                ProblemTest.verifyRelative(testCase,p.exactCurlU(x), ...
                    ProblemTest.curl(@p.exactU,x),'exactCurlU');
                f = ProblemTest.curl(@(z)p.exactCurlU(z)./mu(z),x) + alpha*p.exactU(x);
                ProblemTest.verifyRelative(testCase,p.f(x),f,'f');
                ProblemTest.verifyRelative(testCase,p.curlf(x),ProblemTest.curl(@p.f,x),'curlf');
            else
                ProblemTest.verifyRelative(testCase,p.exactRotU(x), ...
                    ProblemTest.rot(@p.exactU,x),'exactRotU');
                f = ProblemTest.curl(@(z)p.exactRotU(z)./mu(z),x) + alpha*p.exactU(x);
                ProblemTest.verifyRelative(testCase,p.f(x),f,'f');
                ProblemTest.verifyRelative(testCase,p.rotf(x),ProblemTest.rot(@p.f,x),'rotf');
            end
        end

        function testSurfaceExact(testCase,surfaceExact)
            p = surfaceExact();
            x = ProblemTest.surfacePoints(p);
            P = p.normalProjector(x);

            gradU = p.exactGradU(x);
            ProblemTest.verifyRelative(testCase,gradU, ...
                matMul(P,reshape(ProblemTest.jacobian(@p.exactU,x,testCase.h),[],3)), ...
                'tangential gradient');
            w = @(z)matMul(p.normalProjector(z),p.exactGradU(z));
            J = ProblemTest.jacobian(w,x,testCase.h);
            laplaceBeltrami = sum(sum(permute(P,[1 3 2]).*J,2),3);
            ProblemTest.verifyRelative(testCase,p.f(x),-laplaceBeltrami,'f');
        end

        function testSurface(testCase,surface)
            p = surface();
            x = ProblemTest.surfacePoints(p);

            testCase.verifyEqual(p.surfaceFunction(x),zeros(size(x,1),1),'AbsTol',1e-8, ...
                'points on the surface');
            ProblemTest.verifyFD(testCase,p.gradSurfaceFunction(x),@p.surfaceFunction,x, ...
                'gradSurfaceFunction');
            ProblemTest.verifyRelative(testCase,p.d2SurfaceFunction(x), ...
                ProblemTest.jacobian(@p.gradSurfaceFunction,x,testCase.h),'d2SurfaceFunction');
            normal = p.normal(x);
            testCase.verifyEqual(sum(normal.^2,2),ones(size(x,1),1),'AbsTol',testCase.tol);
        end

        function testEigenvalue(testCase,eigenvalue)
            p = eigenvalue();

            val = p.eigenvalue(1:3);

            testCase.verifyTrue(isfinite(val(1)),'first eigenvalue known');
            known = val(isfinite(val));
            testCase.verifyGreaterThanOrEqual(known,0,'nonnegative');
            testCase.verifyTrue(issorted(known),'nondecreasing');
            for k = 1 : 3
                testCase.verifyEqual(p.eigenvalue(k),val(k),sprintf('mode %d',k));
            end
        end

        function testEigenvalueTwist(testCase)
            p = problemScalarEllipticEigenvalueSquareTwist(pi);
            x = ProblemTest.interiorPoints(p.geom,2);
            alpha = p.alpha(x);
            DPhi = p.DPhi(x);

            testCase.verifyEqual(alpha(:,1,2),alpha(:,2,1),'symmetric');
            detAlpha = alpha(:,1,1).*alpha(:,2,2) - alpha(:,1,2).^2;
            testCase.verifyEqual(detAlpha,ones(size(x,1),1),'AbsTol',testCase.tol,'det(alpha)');
            detDPhi = DPhi(:,1,1).*DPhi(:,2,2) - DPhi(:,1,2).*DPhi(:,2,1);
            testCase.verifyEqual(detDPhi,ones(size(x,1),1),'AbsTol',testCase.tol,'det(DPhi)');
            % alpha*DPhi'*DPhi = I
            for i = 1 : 2
                for k = 1 : 2
                    val = 0;
                    for j = 1 : 2
                        val = val + alpha(:,i,j).*sum(DPhi(:,:,j).*DPhi(:,:,k),2);
                    end
                    testCase.verifyEqual(val,(i==k)*ones(size(x,1),1),'AbsTol',testCase.tol, ...
                        sprintf('alpha*DPhi''*DPhi(%d,%d)',i,k));
                end
            end
            ProblemTest.verifyFD(testCase,p.dalpha(x),@p.alpha,x,'dalpha');
            ProblemTest.verifyFD(testCase,DPhi,@p.Phi,x,'DPhi');

            t = linspace(0,1,11)';
            bd = [t,0*t;1+0*t,t;t,1+0*t;0*t,t];
            testCase.verifyEqual(p.Phi(bd),bd,'AbsTol',testCase.tol,'boundary fixed');

            testCase.verifyEqual(p.eigenspace(1),[1 1]);
            testCase.verifyEqual(sortrows(p.eigenspace(2)),[1 2;2 1]);
            testCase.verifyEqual(p.cluster(3),[2 3]);
            m = find(p.eigenvalue(1:40) == 50*pi^2,1);
            testCase.verifyEqual(sortrows(p.eigenspace(m)),[1 7;5 5;7 1]);
            testCase.verifyEqual(p.cluster(m),m:m+2);
            for m = [1 2 4]
                lambda = p.eigenvalue(m);
                for pair = p.eigenspace(m)'
                    u = @(z)p.exactU(z,pair(1),pair(2));
                    gradU = @(z)p.exactGradU(z,pair(1),pair(2));
                    ProblemTest.verifyFD(testCase,gradU(x),u,x,'exactGradU');
                    div = ProblemTest.divergence(@(z)matMul(p.alpha(z),gradU(z)),x);
                    ProblemTest.verifyRelative(testCase,-div,lambda*u(x), ...
                        sprintf('-div(alpha grad u) = lambda u, mode %d',m));
                end
            end
        end

        function testData(testCase,data)
            p = data();
            if isa(p,'ProblemSurfaceTemplate')
                x = ProblemTest.surfacePoints(p);
            else
                x = ProblemTest.interiorPoints(p.geom);
            end

            f = p.f(x);

            testCase.verifySize(f,[size(x,1),size(f,2)]);
            testCase.verifyTrue(all(isfinite(f(:))),'finite right-hand side');
            if ismethod(p,'curlf')
                testCase.verifySize(p.curlf(x),size(f),'curl of the 3d source');
            end
        end

    end

    methods (Static)

        function x = interiorPoints(geom,nrUniform)
            %interiorPoints - Interior points of the elements of the domain
            %   The 2d meshes are refined NRUNIFORM times, once by default.
            if nargin < 2
                nrUniform = 1;
            end
            if size(geom.c4n,2) == 3
                mesh = Mesh3d(geom);
                x = [mesh.point4e([1 1 1]/4);mesh.point4e([0.1 0.2 0.3])];
            else
                mesh = Mesh(geom).refine('uniform',nrUniform);
                x = [mesh.point4e([1 1]/3);mesh.point4e([0.2 0.1])];
            end
        end

        function x = surfacePoints(p)
            %surfacePoints - Centroids of the mesh projected onto the surface
            %   The vertices are avoided since some exact solutions are not
            %   smooth there, e.g. at the poles or reentrant corners.
            mesh = MeshSurface(p.geom,@p.surfaceFunction,@p.gradSurfaceFunction);
            x = projectVertices2Surface(mesh.point4e([1 1]/3),@p.surfaceFunction, ...
                @p.gradSurfaceFunction);
        end

        function J = jacobian(fun,x,h)
            %jacobian - Central differences of the rows of FUN at x
            %   J(:,i,j) is the derivative of the i-th entry of the rows of
            %   FUN in the direction of the j-th coordinate.
            n = size(x,1);
            d = size(x,2);
            J = [];
            for j = 1 : d
                e = zeros(1,d);
                e(j) = h;
                dj = (fun(x + e) - fun(x - e))/(2*h);
                J = cat(3,J,reshape(dj,n,[]));
            end
        end

        function val = divergence(fun,x)
            %divergence - Divergence of the vector field FUN at x
            J = ProblemTest.jacobian(fun,x,ProblemTest.h);
            val = 0;
            for j = 1 : size(x,2)
                val = val + J(:,j,j);
            end
        end

        function val = laplace(gradFun,x)
            %laplace - Laplacian of the components from their gradients
            %   GRADFUN returns the gradients of the components as n x d
            %   (scalar) or n x m x d arrays.
            n = size(x,1);
            d = size(x,2);
            J = ProblemTest.jacobian(@(z)reshape(gradFun(z),size(z,1),[]),x,ProblemTest.h);
            J = reshape(J,n,[],d,d);
            val = 0;
            for j = 1 : d
                val = val + J(:,:,j,j);
            end
        end

        function val = secondDifferences(fun,x,h)
            %secondDifferences - Laplacian of FUN by second differences
            val = 0;
            for j = 1 : size(x,2)
                e = zeros(1,size(x,2));
                e(j) = h;
                val = val + (fun(x + e) - 2*fun(x) + fun(x - e))/h^2;
            end
        end

        function val = curl(fun,x)
            %curl - Curl of the 3d vector field or 2d scalar field FUN at x
            %   In 2d the curl of a scalar field w is (d_2 w, -d_1 w).
            J = ProblemTest.jacobian(fun,x,ProblemTest.h);
            if size(x,2) == 3
                val = [J(:,3,2) - J(:,2,3),J(:,1,3) - J(:,3,1),J(:,2,1) - J(:,1,2)];
            else
                val = [J(:,1,2),-J(:,1,1)];
            end
        end

        function val = rot(fun,x)
            %rot - Rot of the 2d vector or matrix field FUN at x
            %   The rot of a vector field is d_1 v_2 - d_2 v_1 and that of
            %   a matrix field, stored as [11 22 12] if symmetric and as
            %   [11 21 12 22] otherwise, the rots of its rows.
            J = ProblemTest.jacobian(fun,x,ProblemTest.h);
            switch size(J,2)
                case 2
                    val = J(:,2,1) - J(:,1,2);
                case 3
                    val = [J(:,3,1) - J(:,1,2),J(:,2,1) - J(:,3,2)];
                case 4
                    val = [J(:,3,1) - J(:,1,2),J(:,4,1) - J(:,2,2)];
            end
        end

        function verifyRotEpsilon(testCase,p,x)
            %verifyRotEpsilon - Verify the row-wise rot of the strain
            if ismethod(p,'exactRotEpsilon')
                epsilon = @(z)ProblemTest.strain(p.exactGradU(z));
                ProblemTest.verifyRelative(testCase,p.exactRotEpsilon(x), ...
                    ProblemTest.rot(epsilon,x),'exactRotEpsilon');
            end
        end

        function val = strain(gradU)
            %strain - Symmetric gradient [11 22 12] of the 2d gradients
            val = [gradU(:,1,1),gradU(:,2,2),(gradU(:,1,2) + gradU(:,2,1))/2];
        end

        function verifyFD(testCase,val,fun,x,name)
            %verifyFD - Verify derivatives against central differences
            expected = ProblemTest.jacobian(fun,x,ProblemTest.h);
            ProblemTest.verifyRelative(testCase,reshape(val,size(x,1),[]), ...
                reshape(expected,size(x,1),[]),name);
        end

        function verifyRelative(testCase,val,expected,name,relTol)
            %verifyRelative - Verify up to a tolerance relative to the size
            if nargin < 5
                relTol = ProblemTest.fdTol;
            end
            n = size(expected,1);
            expected = reshape(expected,n,[]);
            scale = max(1,max(abs(expected(:))));
            testCase.verifyEqual(reshape(val,n,[]),expected,'AbsTol',relTol*scale,name);
        end

        function val = defines(p,method)
            %defines - True if the class of P defines METHOD itself
            m = findobj(metaclass(p).MethodList,'Name',method);
            val = ~isempty(m) && strcmp(m(1).DefiningClass.Name,class(p));
        end

    end

end
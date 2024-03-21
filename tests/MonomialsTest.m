classdef MonomialsTest < UnitTestCase
    %MonomialsTest - Unit tests of the monomials and their derivatives
    %   Monomials.exponents lists every monomial of total degree at most k
    %   once, nchoosek(k+d,d) rows in lexicographic order. The mass matrix
    %   agrees with Gauss quadrature of the monomial products on the
    %   reference triangle and tetrahedron, and the derivative matrices
    %   applied to the monomials give the falling factorials
    %       d^alpha x^e = e!/(e-alpha)! x^(e-alpha)
    %   (zero for alpha > e) at a point. The constant monomial has
    %   vanishing derivatives, and a scalar coordinate is expanded.
    %
    %   See also Monomials, Pk, Pk3d, PDGk, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        dim = struct('d2',2,'d3',3)
        degree = struct('k0',0,'k1',1,'k2',2,'k4',4,'k6',6)
        order = struct('o1',1,'o2',2,'o3',3,'o4',4)
    end

    methods (Test)

        function testExponents(testCase,dim,degree)
            E = Monomials.exponents(dim,degree);

            testCase.verifySize(E,[nchoosek(degree+dim,dim),dim]);
            testCase.verifyTrue(all(E(:) >= 0) && all(sum(E,2) <= degree),'total degree');
            testCase.verifyEqual(size(unique(E,'rows'),1),size(E,1),'unique');
            testCase.verifyEqual(E,sortrows(E),'lexicographic');
        end

        function testMassMatrix(testCase,dim)
            if dim == 2
                mesh = UnitTestCase.referenceTriangle();
                k = 4;
            else
                mesh = UnitTestCase.referenceTetrahedron();
                k = 3;
            end
            E = Monomials.exponents(dim,k);

            M = Monomials.massMatrix(E);

            f = Monomials.basis(E);
            for m = 1 : size(E,1)
                for n = 1 : size(E,1)
                    val = integrate(@(x,~,~)MonomialsTest.product(f,x,m,n), ...
                        mesh.n4e,1,mesh,2*k);
                    testCase.verifyEqual(M(m,n),val,'AbsTol',testCase.tol, ...
                        sprintf('%s * %s',mat2str(E(m,:)),mat2str(E(n,:))));
                end
            end
        end

        function testDerivatives(testCase,dim,order)
            E = Monomials.exponents(dim,6);
            alpha = Monomials.multiIndices(dim,order);
            x = [0.3 0.7 0.4];
            x = x(1:dim);
            n = size(E,1);

            D = Monomials.derivatives(E,alpha);

            m = MonomialsTest.evaluate(Monomials.basis(E),x).';
            for c = 1 : size(alpha,1)
                rest = E - alpha(c,:);
                expected = prod(factorial(E)./factorial(max(rest,0)).*(rest >= 0),2) ...
                    .*prod(x.^max(rest,0),2);
                testCase.verifyEqual(D((c-1)*n+(1:n),:)*m,expected, ...
                    'AbsTol',testCase.tol*max(1,max(abs(expected))),mat2str(alpha(c,:)));
            end
        end

        function testConstant(testCase,dim)
            E = Monomials.exponents(dim,0);

            testCase.verifyEqual(E,zeros(1,dim));
            testCase.verifyEqual(Monomials.derivatives(E,Monomials.multiIndices(dim,1)), ...
                zeros(dim,1));
            testCase.verifyEqual(Monomials.massMatrix(E),1/factorial(dim),'AbsTol',eps);
        end

        function testScalarCoordinate(testCase)
            f = Monomials.basis(Monomials.exponents(3,3));
            x = [0.1;0.5;0.9];
            y = [0.2;0.4;0.6];

            testCase.verifyEqual(f(x,y,0),f(x,y,zeros(3,1)));
        end

    end

    methods (Static)

        function val = evaluate(f,x)
            %evaluate - Monomials F at the points X, given as rows
            if size(x,2) == 2
                val = f(x(:,1),x(:,2));
            else
                val = f(x(:,1),x(:,2),x(:,3));
            end
        end

        function val = product(f,x,m,n)
            %product - Product of the monomials M and N at the points X
            val = MonomialsTest.evaluate(f,x);
            val = val(:,m).*val(:,n);
        end

    end

end
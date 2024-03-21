classdef QuadratureTest < UnitTestCase
    %QuadratureTest - Unit tests of the Gauss quadrature rules of integrate
    %   integrate(INTEGRAND,PARTS,INDICES,MESH,DEGREE) chooses a Gauss rule
    %   that is exact for polynomials up to degree DEGREE. The tests
    %   integrate all monomials of total degree at most DEGREE over the
    %   unit interval (an edge), the reference triangle, the reference
    %   tetrahedron and the unit square and compare with the exact values
    %       int_0^1 x^a dx              = 1/(a+1)
    %       int_T x^a y^b dx            = a! b! / (a+b+2)!
    %       int_T x^a y^b z^c dx        = a! b! c! / (a+b+c+3)!
    %       int_(0,1)^2 x^a y^b dx      = 1/((a+1)(b+1)).
    %   On the uniformly and locally refined unit square and unit cube the
    %   sum over all elements additionally tests the affine element maps.
    %
    %   See also integrate, UnitTestCase

    % Copyright (c) 2026 by the fem4m authors

    properties (TestParameter)
        degree = struct('d0',0,'d1',1,'d2',2,'d3',3,'d4',4,'d5',5, ...
            'd6',6,'d7',7,'d8',8,'d9',9,'d10',10)
        degree3d = struct('d0',0,'d1',1,'d2',2,'d3',3,'d4',4,'d5',5,'d6',6)
    end

    methods (Test)

        function testInterval(testCase,degree)
            mesh = UnitTestCase.referenceTriangle();
            edge = find(all(mesh.n4ed == [1 2],2));

            for a = 0 : degree
                val = integrate(@(x,~,~)x(:,1).^a,mesh.n4ed(edge,:),edge,mesh,degree);
                testCase.verifyEqual(val,1/(a+1),'AbsTol',testCase.tol, ...
                    sprintf('x^%d',a));
            end
        end

        function testTriangle(testCase,degree)
            mesh = UnitTestCase.referenceTriangle();

            for a = 0 : degree
                for b = 0 : degree-a
                    val = integrate(@(x,~,~)x(:,1).^a.*x(:,2).^b, ...
                        mesh.n4e,1,mesh,degree);
                    exact = factorial(a)*factorial(b)/factorial(a+b+2);
                    testCase.verifyEqual(val,exact,'AbsTol',testCase.tol, ...
                        sprintf('x^%d y^%d',a,b));
                end
            end
        end

        function testTetrahedron(testCase,degree3d)
            mesh = UnitTestCase.referenceTetrahedron();

            for a = 0 : degree3d
                for b = 0 : degree3d-a
                    for c = 0 : degree3d-a-b
                        val = integrate(@(x,~,~)x(:,1).^a.*x(:,2).^b.*x(:,3).^c, ...
                            mesh.n4e,1,mesh,degree3d);
                        exact = factorial(a)*factorial(b)*factorial(c) ...
                            /factorial(a+b+c+3);
                        testCase.verifyEqual(val,exact,'AbsTol',testCase.tol, ...
                            sprintf('x^%d y^%d z^%d',a,b,c));
                    end
                end
            end
        end

        function testSquare(testCase,degree)
            % integrate only needs the coordinates and the area of the
            % element for the tensor Gauss rule on quadrilaterals
            mesh = struct('c4n',[0 0;1 0;1 1;0 1],'area4e',1);

            for a = 0 : degree
                for b = 0 : degree-a
                    val = integrate(@(x,~,~)x(:,1).^a.*x(:,2).^b, ...
                        [1 2 3 4],1,mesh,degree);
                    testCase.verifyEqual(val,1/((a+1)*(b+1)),'AbsTol',testCase.tol, ...
                        sprintf('x^%d y^%d',a,b));
                end
            end
        end

        function testUnitSquareMesh(testCase,degree)
            mesh = PatchTestCase.patchMesh(struct('name','Square'));

            for a = 0 : degree
                b = degree-a;
                val = integrate(@(x,~,~)x(:,1).^a.*x(:,2).^b, ...
                    mesh.n4e,1:mesh.nrElems,mesh,degree);
                testCase.verifyEqual(sum(val),1/((a+1)*(b+1)),'AbsTol',testCase.tol, ...
                    sprintf('x^%d y^%d',a,b));
            end
        end

        function testUnitCubeMesh(testCase,degree3d)
            mesh = UnitTestCase.cubeMesh();

            for a = 0 : degree3d
                for b = 0 : degree3d-a
                    c = degree3d-a-b;
                    val = integrate(@(x,~,~)x(:,1).^a.*x(:,2).^b.*x(:,3).^c, ...
                        mesh.n4e,1:mesh.nrElems,mesh,degree3d);
                    testCase.verifyEqual(sum(val),1/((a+1)*(b+1)*(c+1)), ...
                        'AbsTol',testCase.tol,sprintf('x^%d y^%d z^%d',a,b,c));
                end
            end
        end

    end

end
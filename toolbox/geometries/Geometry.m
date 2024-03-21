classdef Geometry
    %GEOMETRY - Construction of 2d and 3d domain geometries
    %   OBJ = GEOMETRY() creates the coarse triangulation of the unit
    %   square [0,1]^2 with two triangles and Dirichlet boundary
    %   conditions on all edges.
    %
    %   OBJ = GEOMETRY(OPTS) also specifies the options struct OPTS,
    %   whose fields control the construction. Missing fields get their
    %   defaults, and OBJ.opts stores the completed struct.
    % OPTS.name: the domain (default 'Square'). The names of the
    % built-in domains are case-insensitive. Their coarse meshes are
    % built in the centred square (cube) of side 1 and then rotated
    % by OPTS.rotation, scaled by OPTS.length and shifted by
    % OPTS.midpoint; with the defaults this gives:
    %
    % 'Square' ([0,1]^2, up to 4 different boundary conditions)
    %
    %            3
    %         x------x
    %         |      |
    %        4|      |2
    %         |      |
    %         x------x
    %            1
    %
    % 'Lshape' ((-1,1)^2 without its lower right quarter, up to 2
    % different boundary conditions)
    %
    %           2      2
    %        x------x------x
    %        |             |
    %       2|             |2
    %        |             |
    %        x      x------x
    %        |      |   1
    %       2|      |1
    %        |      |
    %        x------x
    %           2
    %
    % 'Slit' ((-1,1)^2 with a slit along [0,1]x{0}, up to 3 different
    % boundary conditions; 1 is the lower and 2 the upper side of the
    % slit)
    %
    %           3      3
    %        x------x------x
    %        |             |
    %       3|             |3
    %        |         2   |
    %        x      x------x
    %        |         1   |
    %       3|             |3
    %        |             |
    %        x------x------x
    %           3      3
    %
    % 'Circle' (unit disk, up to 4 different boundary conditions in
    % the order of 'Square'): the coarse mesh of 'Square' rotated by
    % pi/4 with the corners (1,0), (0,1), (-1,0), (0,-1). Its
    % projectBoundary moves the boundary nodes of a Mesh onto the unit
    % circle after every refinement, and boundaryNormals returns the
    % normals of the circle. The projection and the normals are
    % centred at the origin, so OPTS.midpoint must stay [0,0] and
    % OPTS.length must be a scalar.
    %
    % 'Cube' ([0,1]^3, OPTS.type 'T' only, up to 6 different boundary
    % conditions, in the order of the faces y = min, x = max, y = max,
    % x = min, z = max, z = min)
    %
    % 'Name' (every name and type without a built-in builder) loads
    % the mesh data from the files on the MATLAB path, for example
    % from the subfolders of geometries ('Cooks', 'Fichera',
    % 'Lshape3d', ...); the built-in 'Circle' shadows the files
    % Circle_*.dat, and 'Cube' with OPTS.type 'T' the files Cube_*.dat.
    % OPTS.type, OPTS.reftype and OPTS.bdry are ignored, the boundary
    % sides are given by the files:
    %       Name_c4n.dat
    %       Name_n4e.dat
    %       Name_n4Db.dat (if exists)
    %       Name_n4Nb.dat (if exists)
    %       Name_n4Sb.dat (if exists, side nodes plus fixed axis)
    %       Name_t4e.dat (if exists, 3d bisection labels, see Mesh3d)
    % The error fem4m:Geometry:unknownDomain is raised if there is no
    % builder for OPTS.name and OPTS.type and Name_c4n.dat or
    % Name_n4e.dat is not found, and fem4m:Geometry:invalidN4Sb if
    % Name_n4Sb.dat does not have 3 columns (4 for tetrahedral meshes):
    % the side nodes and the fixed axis.
    %
    % OPTS.type: 'T' (default, triangles, tetrahedra for 'Cube'), 'Q'
    % (quadrilaterals) or 'H' (polygonal mesh of H type for MeshPolyH),
    % case-insensitive; 'Q' and 'H' only for the 2d built-in domains
    % ('Cube' with 'Q' or 'H' loads the tetrahedral mesh of [0,1]^3
    % from the files Cube_*.dat, which the default OPTS.midpoint
    % shifts to [0.5,1.5]^3).
    %
    % OPTS.reftype: pattern of the coarse 'T' meshes of the built-in
    % domains (case-insensitive; other values raise the error
    % fem4m:Geometry:unknownRefType, 'Q' and 'H' ignore it):
    %       'criss' (default) - one diagonal in each square of the 'Q'
    %                           mesh (for 'Square' the diagonal from
    %                           the lower left to the upper right
    %                           corner; for 'Cube' 6 Kuhn tetrahedra)
    %       'cross'           - the other diagonal ('Cube': as 'criss')
    %       'crisscross'      - both diagonals, 4 triangles per square
    %                           ('Cube': 12 tetrahedra around the
    %                           midpoint)
    %       'unionjack'       - 'Square' and 'Circle': 8 triangles,
    %                           the diagonals of the 4 quarter
    %                           squares meet at the midpoint
    %       'corner'          - 'Lshape' and 'Slit' only: all
    %                           diagonals end at the point (0,0)
    %
    % OPTS.bdry: one letter per boundary part, 'D' (Dirichlet), 'N'
    % (Neumann) or 'S' (symmetry, the displacement component normal to
    % the side is fixed), e.g. 'DSNSDS' for the six faces of 'Cube'. A
    % single letter applies to all parts (default 'D'). The fixed axis
    % of an 'S' side is the coordinate axis parallel to its normal, or
    % 0 for a general normal; a nonzero OPTS.rotation sets it to 0 for
    % all 2d sides (e.g. for 'Circle').
    %
    % OPTS.length: scaling factor of the coordinates, a scalar, [LX LY]
    % or [LX LY LZ] (default sqrt(2) for 'Circle', 2 for 'Lshape' and
    % 'Slit', 1 otherwise).
    %
    % OPTS.midpoint: shift added to the coordinates after the scaling,
    % a scalar, 1-by-2 or 1-by-3 (default [0.5,0.5] for 'Square',
    % [0.5,0.5,0.5] for 'Cube', [0,0] otherwise).
    %
    % OPTS.rotation: angle of the rotation about the origin (in 3d
    % about the z axis) applied before the scaling (default pi/4 for
    % 'Circle', 0 otherwise).
    %
    % OPTS.slitdelta: opening of the slit of 'Slit' (default 0). For
    % 'T' and 'Q' the right end of the lower slit side lies
    % OPTS.slitdelta*LY below that of the upper side (LY = 2 by
    % default); 'H' moves the nodes of the lower slit side down by
    % half of that.
    %
    % OPTS.projectBoundary: function handle MESH = F(MESH) stored in
    % projectBoundary, which Mesh calls after every refinement and
    % coarsening (default: project2Circle with the radius
    % sqrt(2)*OPTS.length/2 for 'Circle', [] otherwise).
    %
    %   GEOMETRY functions:
    %       squareT        - Build coarse triangular mesh of the unit square
    %       squareQ        - Build coarse quadrilateral mesh of the unit square
    %       squareH        - Build coarse H-type mesh of the unit square
    %       lshapeT        - Build coarse triangular mesh of the L-shape
    %       lshapeQ        - Build coarse quadrilateral mesh of the L-shape
    %       lshapeH        - Build coarse H-type mesh of the L-shape
    %       slitT          - Build coarse triangular mesh of the slit domain
    %       slitQ          - Build coarse quadrilateral mesh of the slit domain
    %       slitH          - Build coarse H-type mesh of the slit domain
    %       cubeT          - Build coarse tetrahedral mesh of the unit cube
    %       project2Circle - Project boundary nodes onto a circle
    %
    %   GEOMETRY properties:
    %       c4n             - Coordinates of nodes
    %       n4e             - Node ids for elements (cell array for 'H')
    %       n4Db            - Node pairs (3d: triples) of the Dirichlet
    %                         boundary sides
    %       n4Nb            - Node pairs (3d: triples) of the Neumann
    %                         boundary sides
    %       n4Sb            - Node pairs (3d: triples) of the symmetry
    %                         boundary sides plus the fixed axis
    %                         (2d: 0/1/2, 3d: 0/1/2/3)
    %       t4e             - Bisection labels of 3d meshes, [] for the
    %                         longest edges
    %       n4Bdry          - Cell array of the node pairs (3d: triples)
    %                         of the sides per boundary part
    %       axis4Bdry       - Cell array of the fixed axes (1, 2 or 0)
    %                         of the sides per boundary part before
    %                         the rotation, [] for 'Cube' and for
    %                         meshes loaded from files
    %       projectBoundary - Function handle MESH = F(MESH) moving the
    %                         boundary nodes onto a curved boundary,
    %                         [] for none
    %       boundaryNormals - Function handle of the unit outer normals
    %                         at boundary points (N-by-2), [] except
    %                         for 'Circle'
    %       opts            - Geometry options used for construction
    %
    %   See also Mesh, Mesh3d, MeshSurface, MeshPoly, MeshPolyH

    % Copyright (c) 2026 by the fem4m authors

    properties
        c4n             % coordinates of nodes
        n4e             % node ids for elements (cell array for 'H')
        n4Db            % node pairs (3d: triples) of the Dirichlet sides
        n4Nb            % node pairs (3d: triples) of the Neumann sides
        n4Sb            % node pairs (3d: triples) of symmetry sides, fixed axis
        t4e             % bisection labels of 3d meshes, [] for longest edges
        n4Bdry          % cell array of the side nodes per boundary part
        axis4Bdry       % cell array of the fixed axes per boundary part
        projectBoundary % function handle MESH = F(MESH), [] for none
        boundaryNormals % function handle of the unit outer normals
    end

    properties(GetAccess=public,SetAccess=private)
        opts            % geometry options used for construction
    end

    methods
        function obj = Geometry(opts)
            %GEOMETRY - Create a geometry object from options
            if nargin < 1
                opts = [];
            end
            if ~isfield(opts,'name')
                opts.name = 'square';
            end
            if ~isfield(opts,'type')
                opts.type = 'T';
            end
            if ~isfield(opts,'reftype')
                opts.reftype = 'criss';
            end
            if ~isfield(opts,'midpoint')
                if strcmpi(opts.name,'square')
                    opts.midpoint = [0.5,0.5];
                elseif strcmpi(opts.name,'cube')
                    opts.midpoint = [0.5,0.5,0.5];
                else
                    opts.midpoint = [0,0];
                end
            end
            if ~isfield(opts,'bdry')
                opts.bdry = 'D';
            end
            if ~isfield(opts,'length')
                if strcmpi(opts.name,'square')
                    opts.length = 1;
                elseif strcmpi(opts.name,'circle')
                    opts.length = sqrt(2);
                elseif strcmpi(opts.name,'lshape')
                    opts.length = 2;
                elseif strcmpi(opts.name,'slit')
                    opts.length = 2;
                else
                    opts.length = 1;
                end
            end
            if ~isfield(opts,'rotation')
                if strcmpi(opts.name,'circle')
                    opts.rotation = pi/4;
                else
                    opts.rotation = 0;
                end
            end
            if ~isfield(opts,'slitdelta')
                opts.slitdelta = 0;
            end
            obj.opts = opts;

            if isfield(opts,'projectBoundary')
                obj.projectBoundary = opts.projectBoundary;
            elseif strcmpi(opts.name,'circle')
                obj.projectBoundary = @(mesh)obj.project2Circle(mesh,sqrt(2*opts.length^2)/2);
            else
                obj.projectBoundary = [];
            end

            if strcmpi(opts.name,'circle')
                fname = 'square';
                obj.boundaryNormals = @(z)(z./vecnorm(z,2,2));
            else
                fname = lower(opts.name);
            end

            if ismethod(obj,[fname,upper(opts.type)])

                obj = obj.([fname,upper(opts.type)])(opts);

                if isscalar(opts.bdry)
                    opts.bdry = repmat(opts.bdry,numel(obj.n4Bdry),1);
                end
                for curEd = 1 : numel(obj.n4Bdry)
                    n4xb = ['n4',upper(opts.bdry(curEd)),'b'];
                    obj.(n4xb) = [obj.(n4xb);obj.n4Bdry{curEd}];
                end

                % symmetry boundary: the fixed axis in the last column
                % (2d: 1 or 2, 3d: 1, 2 or 3, 0 for a general normal)
                if ~isempty(obj.n4Sb)
                    axis4Sb = zeros(0,1);
                    for curEd = 1 : numel(obj.n4Bdry)
                        if upper(opts.bdry(curEd)) == 'S'
                            if numel(obj.axis4Bdry) >= curEd && ~isempty(obj.axis4Bdry{curEd})
                                axis4Sb = [axis4Sb;obj.axis4Bdry{curEd}(:)]; %#ok<AGROW>
                            else
                                axis4Sb = [axis4Sb;fixedAxis(obj.c4n,obj.n4Bdry{curEd})]; %#ok<AGROW>
                            end
                        end
                    end
                    obj.n4Sb(:,end+1) = axis4Sb;
                end

            elseif exist([opts.name,'_c4n.dat'],'file') && exist([opts.name,'_n4e.dat'],'file')

                obj.c4n   = load([opts.name,'_c4n.dat']);
                obj.n4e   = load([opts.name,'_n4e.dat']);
                if exist([opts.name,'_n4Db.dat'],'file')
                    obj.n4Db = load([opts.name,'_n4Db.dat']);
                else
                    obj.n4Db = [];
                end
                if exist([opts.name,'_n4Nb.dat'],'file')
                    obj.n4Nb  = load([opts.name,'_n4Nb.dat']);
                else
                    obj.n4Nb = [];
                end
                if exist([opts.name,'_n4Sb.dat'],'file')
                    obj.n4Sb  = load([opts.name,'_n4Sb.dat']);
                else
                    obj.n4Sb = [];
                end
                % symmetry sides (node triples of tetrahedra, otherwise
                % node pairs) with the fixed axis in the last column
                nrNodes4Sb = 2 + (size(obj.c4n,2) == 3 && size(obj.n4e,2) == 4);
                if ~isempty(obj.n4Sb) && size(obj.n4Sb,2) ~= nrNodes4Sb+1
                    error('fem4m:Geometry:invalidN4Sb', ...
                        '%s_n4Sb.dat must have %d columns: the side nodes and the fixed axis.', ...
                        opts.name,nrNodes4Sb+1);
                end
                if exist([opts.name,'_t4e.dat'],'file')
                    obj.t4e = load([opts.name,'_t4e.dat']);
                end

            else
                error('fem4m:Geometry:unknownDomain', ...
                    'Create geometry: domain ''%s'' for %s-meshes not implemented and files %s_c4n.dat, %s_n4e.dat not found!', ...
                    lower(opts.name),upper(opts.type),opts.name,opts.name);
            end

            %rotation
            if opts.rotation ~= 0
                % rotation about the z axis in 3d
                [theta,rho] = cart2pol(obj.c4n(:,1),obj.c4n(:,2));
                [x,y] = pol2cart(theta+opts.rotation,rho);
                obj.c4n(:,1:2) = [x,y];
                if ~isempty(obj.n4Sb) && size(obj.c4n,2) == 3 && size(obj.n4e,2) == 4
                    obj.n4Sb(:,4) = fixedAxis(obj.c4n,obj.n4Sb(:,1:3));
                elseif ~isempty(obj.n4Sb)
                    obj.n4Sb(:,3) = 0;
                end
            end

            % linear transform
            if isempty(opts.length)
                opts.length = 1;
            end
            if isscalar(opts.length)
                lengthX = opts.length;
                lengthY = opts.length;
                lengthZ = opts.length;
            elseif numel(opts.length)==2
                lengthX = opts.length(1);
                lengthY = opts.length(2);
                lengthZ = 1;
            else
                lengthX = opts.length(1);
                lengthY = opts.length(2);
                lengthZ = opts.length(3);
            end
            if isempty(opts.midpoint)
                opts.midpoint = 0;
            end
            if isscalar(opts.midpoint)
                midX = opts.midpoint;
                midY = opts.midpoint;
                midZ = opts.midpoint;
            elseif numel(opts.midpoint)==2
                midX = opts.midpoint(1);
                midY = opts.midpoint(2);
                midZ = 0;
            else
                midX = opts.midpoint(1);
                midY = opts.midpoint(2);
                midZ = opts.midpoint(3);
            end

            obj.c4n(:,1) = obj.c4n(:,1)*lengthX + midX;
            obj.c4n(:,2) = obj.c4n(:,2)*lengthY + midY;
            if size(obj.c4n,2)==3
                obj.c4n(:,3) = obj.c4n(:,3)*lengthZ + midZ;
            end

        end

        %% Square

        function obj = squareT(obj,opts)
            %squareT - Build coarse triangular mesh of the unit square
            %   OBJ = squareT(OBJ,OPTS) sets the node coordinates and
            %   elements of a coarse triangulation of the square
            %   [-1/2,1/2]^2 according to
            %   OPTS.reftype and lists its four boundary parts with
            %   their fixed axes in n4Bdry and axis4Bdry. The
            %   constructor calls it for 'Square' and 'Circle' with
            %   OPTS.type 'T'.
            obj.c4n = [-0.5 -0.5; 0.5 -0.5; 0.5 0.5; -0.5 0.5];
            if strcmpi(opts.reftype,'criss')
                obj.n4e = [3 1 2; 1 3 4];
            elseif strcmpi(opts.reftype,'cross')
                obj.n4e = [2 4 1; 4 2 3];
            elseif strcmpi(opts.reftype,'crisscross')
                obj.n4e = [1 2 5; 2 3 5; 3 4 5; 4 1 5];
                obj.c4n = [obj.c4n;0,0];
            elseif strcmpi(opts.reftype,'unionjack')
                obj.n4e = [5 1 2;3 5 2;5 3 6;9 5 6;5 9 8;7 5 8;5 7 4;1 5 4];
                obj.c4n = [0 0;0.5 0;1 0;0 0.5;0.5 0.5;1 0.5;0 1;0.5 1;1 1] - 0.5;
            else
                error('fem4m:Geometry:unknownRefType','reftype ''%s'' not implemented !!!',opts.reftype);
            end
            if strcmpi(opts.reftype,'unionjack')
                obj.n4Bdry = { [1 2;2 3] ; [3 6;6 9]; [9 8;8 7]; [7 4;4 1] };
                obj.axis4Bdry = {[2;2] ; [1;1] ; [2;2] ; [1;1] };
            else
                obj.n4Bdry = { [1 2]; [2 3]; [3 4]; [4 1]};
                obj.axis4Bdry = { 2 ; 1 ; 2 ; 1};
            end
        end

        function obj = squareQ(obj,~)
            %squareQ - Build coarse quadrilateral mesh of the unit square
            %   OBJ = squareQ(OBJ,OPTS) sets a single quadrilateral
            %   element for the square [-1/2,1/2]^2 and lists its four
            %   boundary parts with their fixed axes; OPTS is ignored.
            obj.c4n = [-0.5 -0.5; 0.5 -0.5; 0.5 0.5; -0.5 0.5];
            obj.n4e = [1 2 3 4];
            obj.n4Bdry = {[1 2]; [2 3]; [3 4]; [4 1]};
            obj.axis4Bdry = { 2 ; 1 ; 2 ; 1};
        end

        function obj = squareH(obj,~)
            %squareH - Build coarse H-type mesh of the unit square
            %   OBJ = squareH(OBJ,OPTS) sets a coarse polygonal mesh of
            %   type H (one hexagon, four pentagons and two
            %   quadrilaterals, n4e a cell array) for the square
            %   [-1/2,1/2]^2 and lists its four boundary parts with
            %   their fixed axes; OPTS is ignored.
            obj.c4n = [0 0; 0 0.5; 0 1;
                1/6 0; 1/3 1/4; 1/6 1/2; 1/3 3/4; 1/6 1;
                2/3 1/4; 2/3 3/4; 5/6 0; 5/6 1/2; 5/6 1;
                1 0; 1 1/2; 1 1]-[1/2 1/2];
            obj.n4e = { [     5     9    12    10     7     6];
                [     1     4     5     6     2];
                [    11    14    15    12     9];
                [     2     6     7     8     3];
                [    12    15    16    13    10];
                [     8     7    10    13];
                [    11     9     5     4]};
            obj.n4Bdry = {[1 4; 4 11; 11 14];
                [14 15; 15 16];
                [16 13; 13 8; 8 3];
                [3 2; 2 1]};
            obj.axis4Bdry = { [2;2;2] ; [1;1] ; [2;2;2] ; [1;1]};
        end

        %% Lshape

        function obj = lshapeT(obj,opts)
            %lshapeT - Build coarse triangular mesh of the L-shape
            %   OBJ = lshapeT(OBJ,OPTS) sets the node coordinates and
            %   elements of a coarse triangulation of the L-shape
            %   (-1/2,1/2)^2 without its lower right quarter according
            %   to OPTS.reftype and lists its two boundary parts with
            %   their fixed axes.
            obj.c4n = [0 -1/2; 0 0; 1/2 0; 1/2 1/2; 0 1/2;
                -1/2 1/2; -1/2 0; -1/2 -1/2];
            if strcmpi(opts.reftype,'criss')
                obj.n4e = [2 8 1; 8 2 7; 5 7 2; 7 5 6; 4 2 3; 2 4 5];
            elseif strcmpi(opts.reftype,'cross')
                obj.n4e = [1 7 8; 7 1 2; 2 6 7; 6 2 5; 3 5 2; 5 3 4];
            elseif strcmpi(opts.reftype,'crisscross')
                obj.n4e = [1 2 9; 2 7 9; 7 8 9; 8 1 9; 2 5 10; 5 6 10; 6 7 10;
                    7 2 10; 3 4 11; 4 5 11; 5 2 11; 2 3 11];
                obj.c4n = [obj.c4n;-1/4 -1/4; -1/4 1/4; 1/4 1/4];
            elseif strcmpi(opts.reftype,'corner')
                obj.n4e = [2 8 1; 8 2 7; 2 6 7; 6 2 5; 4 2 3; 2 4 5];
            else
                error('fem4m:Geometry:unknownRefType','reftype ''%s'' not implemented !!!',opts.reftype);
            end
            obj.n4Bdry = {[1 2; 2 3];
                [3 4; 4 5; 5 6; 6 7; 7 8; 8 1]};
            obj.axis4Bdry = { [1; 2]; [1; 2; 2; 1; 1; 2]};
        end

        function obj = lshapeQ(obj,~)
            %lshapeQ - Build coarse quadrilateral mesh of the L-shape
            %   OBJ = lshapeQ(OBJ,OPTS) sets three square elements for
            %   the L-shape (-1/2,1/2)^2 without its lower right quarter
            %   and lists its two boundary parts with their fixed axes;
            %   OPTS is ignored.
            obj.c4n = [0 -1/2; 0 0; 1/2 0; 1/2 1/2; 0 1/2;
                -1/2 1/2; -1/2 0; -1/2 -1/2];
            obj.n4e = [8 1 2 7; 7 2 5 6; 2 3 4 5];
            obj.n4Bdry = {[1 2; 2 3];
                [3 4; 4 5; 5 6; 6 7; 7 8; 8 1]};
            obj.axis4Bdry = { [1; 2]; [1; 2; 2; 1; 1; 2]};
        end

        function obj = lshapeH(obj,~)
            %lshapeH - Build coarse H-type mesh of the L-shape
            %   OBJ = lshapeH(OBJ,OPTS) merges three once refined H-type
            %   meshes of squares (MeshPolyH) into the L-shape
            %   (-1/2,1/2)^2 without its lower right quarter and lists
            %   its two boundary parts with their fixed axes; OPTS is
            %   ignored.
            geom.name = 'square';
            geom.type = 'H';
            geom1 = Geometry(geom);
            mesh1 = MeshPolyH(geom1.c4n,geom1.n4e,geom1.n4Db,[]);
            mesh1 = mesh1.refine(ones(mesh1.nrElems,1));
            geom.midpoint = [-0.5,0.5];
            geom2 = Geometry(geom);
            mesh2 = MeshPolyH(geom2.c4n,geom2.n4e,geom2.n4Db,[]);
            mesh2 = mesh2.refine(ones(mesh2.nrElems,1));
            geom.midpoint = [-0.5,-0.5];
            geom3 = Geometry(geom);
            mesh3 = MeshPolyH(geom3.c4n,geom3.n4e,geom3.n4Db,[]);
            mesh3 = mesh3.refine(ones(mesh3.nrElems,1));
            mesh = mesh1+mesh2+mesh3;

            c4nH = mesh.c4n;
            n4eH = mesh.n4e;
            n4DbH = mesh.n4Db;

            n4eH{25} = [64,66,68,62,103,105];
            n4eH{32} = [71,109,125,51,73];
            n4eH{19} = [7,60,76,41];
            n4eH{7} = [33,37,41,76,80,84];
            n4eH{10} = [39,5,33,84,59,74];
            n4eH{5} = [31,35,39,74,78,82];
            n4eH{9} = [114,1,3,31,82,58,118];
            n4eH([47,65,34,31,43,29,41,57])=[];
            n4DbH(15,:) = [7,60];
            n4DbH(23,:) = [71,109];
            n4DbH([24,31],:) = [];

            [c4nH,n4eH,n4DbH,~] = MeshPoly.removeGhostNodes(c4nH,n4eH,n4DbH,[]);

            obj.c4n = 0.5*c4nH;
            obj.n4e = n4eH;
            indx1 = [28 31 29 32 35 12 1 3 5 14];
            indx2 = setdiff(1:size(n4DbH,1),indx1);
            obj.n4Bdry = { n4DbH(indx1,:) ; n4DbH(indx2,:) };
            axis4DbH = fixedAxis(obj.c4n,n4DbH);
            obj.axis4Bdry = { axis4DbH(indx1) ; axis4DbH(indx2) };
        end

        %% Slit

        function obj = slitT(obj,opts)
            %slitT - Build coarse triangular mesh of the slit domain
            %   OBJ = slitT(OBJ,OPTS) sets the node coordinates and
            %   elements of a coarse triangulation of the square
            %   (-1/2,1/2)^2 with a slit from (0,0) to (1/2,0) according
            %   to OPTS.reftype and lists its three boundary parts with
            %   their fixed axes. The right end of the lower slit side
            %   is (1/2,-OPTS.slitdelta).
            delta = opts.slitdelta;
            obj.c4n = [-1/2 -1/2; 1/2 -1/2; -1/2 1/2; 1/2 1/2; 0 -1/2;
                -1/2 0; 0 0; 1/2 -delta; 0 1/2; 1/2 0];
            if strcmpi(opts.reftype,'criss')
                obj.n4e = [4 7 10; 7 4 9; 9 6 7; 6 9 3;
                    7 1 5; 1 7 6; 8 5 2; 5 8 7];
            elseif strcmpi(opts.reftype,'cross')
                obj.n4e = [6 5 7; 10 9 7; 2 7 5; 3 7 9;
                    7 3 6; 7  2 8; 5 6 1; 9 10 4];
            elseif strcmpi(opts.reftype,'crisscross')
                obj.n4e = [1 5 11; 5 7 11; 7 6 11; 6 1 11;
                    6 7 12; 7 9 12; 9 3 12; 3 6 12;
                    7 10 13; 10 4 13; 4 9 13; 9 7 13;
                    5 2 14; 2 8 14; 8 7 14; 7 5 14];
                obj.c4n = [obj.c4n;-1/4 -1/4; -1/4 1/4; 1/4 1/4; 1/4 -1/4];
            elseif strcmpi(opts.reftype,'corner')
                obj.n4e = [7 1 5; 1 7 6; 3 7 9; 7 3 6;
                    4 7 10; 7 4 9; 7 2 8; 2 7 5];
            else
                error('fem4m:Geometry:unknownRefType','reftype ''%s'' not implemented !!!',opts.reftype);
            end
            obj.n4Bdry = {[8 7]; [7 10];
                [10 4; 4 9; 9 3; 3 6; 6 1; 1 5; 5 2; 2 8]};
            obj.axis4Bdry = { 2; 2; [1;2;2;1;1;2;2;1]};

        end

        function obj = slitQ(obj,opts)
            %slitQ - Build coarse quadrilateral mesh of the slit domain
            %   OBJ = slitQ(OBJ,OPTS) sets four quadrilateral elements
            %   for the square (-1/2,1/2)^2 with a slit from (0,0) to
            %   (1/2,0) and lists its three boundary parts with their
            %   fixed axes. The right end of the lower slit side is
            %   (1/2,-OPTS.slitdelta).
            delta = opts.slitdelta;
            obj.c4n = [-1/2 -1/2; 1/2 -1/2; -1/2 1/2; 1/2 1/2; 0 -1/2;
                -1/2 0; 0 0; 1/2 -delta; 0 1/2; 1/2 0];
            obj.n4e = [1 5 7 6; 6 7 9 3; 5 2 8 7; 7 10 4 9];
            obj.n4Bdry = {[8 7]; [7 10];
                [10 4; 4 9; 9 3; 3 6; 6 1; 1 5; 5 2; 2 8]};
            obj.axis4Bdry = { 2; 2; [1;2;2;1;1;2;2;1]};
        end

        function obj = slitH(obj,opts)
            %slitH - Build coarse H-type mesh of the slit domain
            %   OBJ = slitH(OBJ,OPTS) merges four once refined H-type
            %   meshes of squares (MeshPolyH) into the square
            %   (-1/2,1/2)^2 with a slit from (0,0) to (1/2,0) and lists
            %   its three boundary parts with their fixed axes. The
            %   nodes of the lower slit side are moved down by
            %   OPTS.slitdelta/2.
            delta = 1e-6;
            geom.name = 'square';
            geom.type = 'H';
            geom1 = Geometry(geom);
            mesh1 = MeshPolyH(geom1.c4n,geom1.n4e,geom1.n4Db,[]);
            mesh1 = mesh1.refine(ones(mesh1.nrElems,1));
            geom.midpoint = [-0.5,0.5];
            geom2 = Geometry(geom);
            mesh2 = MeshPolyH(geom2.c4n,geom2.n4e,geom2.n4Db,[]);
            mesh2 = mesh2.refine(ones(mesh2.nrElems,1));
            geom.midpoint = [-0.5,-0.5];
            geom3 = Geometry(geom);
            mesh3 = MeshPolyH(geom3.c4n,geom3.n4e,geom3.n4Db,[]);
            mesh3 = mesh3.refine(ones(mesh3.nrElems,1));
            geom.midpoint = [0.5,-0.5];
            geom4 = Geometry(geom);
            mesh4 = MeshPolyH(geom4.c4n,geom4.n4e,geom4.n4Db,[]);
            mesh4 = mesh4.refine(ones(mesh4.nrElems,1));
            indx = [7 17 15 12 14];
            mesh4.c4n(indx,2) = mesh4.c4n(indx,2)-delta;
            mesh = mesh1+mesh2+mesh3+mesh4;

            c4nH = mesh.c4n;
            n4eH = mesh.n4e;
            n4DbH = mesh.n4Db;

            n4eH{25} = [64,66,68,62,103,105];
            n4eH{32} = [71,109,125,51,73];
            n4eH{19} = [7,60,76,41];
            n4eH{7} = [33,37,41,76,80,84];
            n4eH{10} = [39,5,33,84,59,74];
            n4eH{5} = [31,35,39,74,78,82];
            n4eH{9} = [167,137,1,3,31,82,58,118];
            n4eH{76} = [126,159,163,167,118,122];
            n4eH{79} = [116,165,135,159,126,99];
            n4eH{74} = [124,157,161,165,116,120];
            n4eH{78} = [98,133,157,124];
            n4eH([47,65,34,31,43,29,41,57,88,54,66,52,64])=[];

            n4DbH(15,:) = [7,60];
            n4DbH(23,:) = [71,109];
            n4DbH(35,:) = [98,133];
            n4DbH([24,29,47],:) = [];

            [c4nH,n4eH,n4DbH,~] = MeshPoly.removeGhostNodes(c4nH,n4eH,n4DbH,[]);

            indx = [126,136,134,131,133];
            c4nH(indx,2) = c4nH(indx,2)+delta;
            % fixed axes of the slit of height 0, as in slitT and slitQ
            axis4DbH = fixedAxis(c4nH,n4DbH);
            c4nH(indx,2) = c4nH(indx,2)-opts.slitdelta;
            obj.c4n = 0.5*c4nH;
            obj.n4e = n4eH;
            indx1 = [45, 35, 37, 39, 47];
            indx2 = [12, 1, 3, 5, 14];
            indx3 = setdiff(1:size(n4DbH,1),[indx1,indx2]);
            obj.n4Bdry = { n4DbH(indx1,:); n4DbH(indx2,:); n4DbH(indx3,:)};
            obj.axis4Bdry = { axis4DbH(indx1); axis4DbH(indx2); axis4DbH(indx3)};
        end

        %% cube
        function obj = cubeT(obj,opts)
            %cubeT - Build coarse tetrahedral mesh of the unit cube
            %   OBJ = cubeT(OBJ,OPTS) sets the node coordinates,
            %   elements and bisection labels t4e of a coarse
            %   tetrahedral mesh of the cube [-1/2,1/2]^3 according to
            %   OPTS.reftype ('criss' or 'cross': 6 Kuhn tetrahedra,
            %   'crisscross': 12 tetrahedra around the midpoint) and
            %   lists its six boundary parts (one per face). It sets no
            %   axis4Bdry; the constructor computes the fixed axes of
            %   symmetry faces from their normals.
            if strcmpi(opts.reftype,'criss') || strcmpi(opts.reftype,'cross')
                obj.c4n =    [
                    -0.5   -0.5   -0.5
                    -0.5    0.5   -0.5
                    0.5   -0.5   -0.5
                    0.5    0.5   -0.5
                    -0.5   -0.5    0.5
                    -0.5    0.5    0.5
                    0.5   -0.5    0.5
                    0.5    0.5    0.5];
                % bisection order: the refinement edge n4e(:,[1 4]) is
                % the diagonal 4-5 (Kuhn simplices), type 0
                obj.n4e = [5 1 3 4
                    4 2 1 5
                    5 6 2 4
                    4 8 6 5
                    4 3 7 5
                    5 7 8 4];
                obj.t4e = zeros(6,1);
                obj.n4Bdry = { [ 3 5 1; 5 3 7 ];
                    [ 4 7 3; 7 4 8 ];
                    [ 4 6 8; 6 4 2];
                    [ 2 5 6; 5 2 1];
                    [ 5 8 6; 8 5 7];
                    [ 1 4 3; 4 1 2]};
            elseif strcmpi(opts.reftype,'crisscross')
                obj.c4n =    [-0.5   -0.5   -0.5
                    -0.5    0.5   -0.5
                    0.5   -0.5   -0.5
                    0.5    0.5   -0.5
                    -0.5   -0.5    0.5
                    -0.5    0.5    0.5
                    0.5   -0.5    0.5
                    0.5    0.5    0.5
                    0         0         0];
                % bisection order: the refinement edge n4e(:,[1 4]) is a
                % face diagonal, type 1
                obj.n4e = [4 9 3 1
                    1 9 2 4
                    1 9 3 7
                    7 9 5 1
                    7 9 3 4
                    4 9 8 7
                    2 9 1 5
                    5 9 6 2
                    8 9 4 2
                    2 9 6 8
                    8 9 6 5
                    5 9 7 8];
                obj.t4e = ones(12,1);
                % one part per cube face, in the order of the criss cube
                obj.n4Bdry = { [ 1 3 7; 1 7 5];
                    [ 3 4 7; 4 8 7];
                    [ 2 8 4; 2 6 8];
                    [ 1 5 2; 2 5 6];
                    [ 5 8 6; 5 7 8];
                    [ 1 2 4; 1 4 3]};
            else
                error('fem4m:Geometry:unknownRefType','reftype ''%s'' not implemented !!!',opts.reftype);
            end
        end

    end

    methods (Static)

        function mesh = project2Circle(mesh,r)
            %project2Circle - Project boundary nodes onto a circle
            %   MESH = project2Circle(MESH,R) scales the coordinates of
            %   the nodes of all boundary edges (n4Db, n4Nb and n4Sb) of
            %   the 2d mesh MESH radially onto the circle with radius R
            %   centred at the origin.
            n4Bd = [mesh.n4Db;mesh.n4Nb];
            if ~isempty(mesh.n4Sb)
                n4Bd = [n4Bd;mesh.n4Sb(:,1:2)];
            end
            for curEdge = 1 : size(n4Bd,1)
                for k = 1 : 2
                    x = mesh.c4n(n4Bd(curEdge,k),1);
                    y = mesh.c4n(n4Bd(curEdge,k),2);
                    mesh.c4n(n4Bd(curEdge,k),:) = mesh.c4n(n4Bd(curEdge,k),:).*(r/sqrt(x^2+y^2));
                end
            end
        end

    end
end

function axis4s = fixedAxis(c4n,n4s)
%fixedAxis - Fixed axis of symmetry boundary sides
%   AXIS4S = fixedAxis(C4N,N4S) returns for the boundary sides with the
%   node ids N4S (node pairs of edges in 2d, node triples of faces in 3d)
%   the number k of the fixed displacement component if the side normal
%   is parallel to the k-th unit vector, and 0 (general normal)
%   otherwise. Normal components off the axis may be up to sqrt(eps)
%   times the length of the normal, which allows the rounding errors of
%   transformed or far-off coordinates.
c1 = c4n(n4s(:,1),:);
if size(n4s,2) == 2
    t = c4n(n4s(:,2),:) - c1;
    normal = [t(:,2),-t(:,1)];
else
    normal = cross(c4n(n4s(:,2),:)-c1,c4n(n4s(:,3),:)-c1,2);
end
small = abs(normal) <= sqrt(eps)*vecnorm(normal,2,2);
axis4s = ~small*(1:size(normal,2))';
axis4s(sum(small,2) ~= size(normal,2)-1) = 0;
end
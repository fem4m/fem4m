classdef fem4mResultTest < UnitTestCase
    %fem4mResultTest - File names and contents of saved level-wise results
    %   The level data are the numbers of degrees of freedom and the H1
    %   and L2 errors of P2 on three uniformly refined meshes of the
    %   Poisson problem of problemScalarEllipticSquareOscExact. The tests
    %   check the file names that plotApp splits into eight slots, the
    %   p-encoding of numbers and the texts plotApp shows for them, the
    %   merge of the labels, the contents of the saved struct p for
    %   levels given as name-value arguments or as struct, and the
    %   results files of the app modules.
    %
    %   See also fem4mResult, fem4mModuleTemplate, ConvergenceTest

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = private)
        Levels  % Struct of nrDoF, H1Error and L2Error per level of P2
        Folder  % Temporary folder of the saved files
    end

    methods (TestClassSetup)

        function computeLevels(testCase)
            p = problemScalarEllipticSquareOscExact();
            mesh = Mesh(p.geom);
            method = struct('name','Pk','degree',2);
            levels = struct('nrDoF',zeros(1,3),'H1Error',zeros(1,3),'L2Error',zeros(1,3));
            for level = 1 : 3
                mesh = mesh.refine('uniform');
                levels.nrDoF(level) = Pk(2,mesh).nrDoF;
                err = ConvergenceTest.solve(p,mesh,method);
                levels.H1Error(level) = err(1);
                levels.L2Error(level) = err(2);
            end
            testCase.Levels = levels;
        end

    end

    methods (TestMethodSetup)

        function temporaryFolder(testCase)
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end

    end

    methods (Test)

        function testFileName(testCase)
            testCase.verifyEqual(fem4mResultTest.pkResult().fileName(), ...
                'Square_Pk_uniform_residual_direct_2_0_10000');
        end

        function testEmptySlots(testCase)
            testCase.verifyError(@()fem4mResult().fileName(), ...
                'fem4m:fem4mResult:emptySlot');
            for j = 1 : numel(fem4mResult.Slots)
                r = fem4mResultTest.pkResult();
                r.(fem4mResult.Slots{j}) = [];
                testCase.verifyError(@()r.fileName(),'fem4m:fem4mResult:emptySlot', ...
                    fem4mResult.Slots{j});
            end
        end

        function testFileNameSlots(testCase)
            r = fem4mResult(problem='Square',fe='Pk',refine='adaptive', ...
                estimator="none",solver='pcg',param1=2,param2=0.3, ...
                maxNrDoF=150000);
            name = r.fileName();
            testCase.verifyEqual(name, ...
                'Square_Pk_adaptive_none_pcg_2_0p3_150000');
            % as split by plotApp and mat2data
            parts = strsplit(name,'_');
            testCase.verifyNumElements(parts,8);
            testCase.verifyEqual(fem4mResult.decodeNumber(parts{7}),0.3);
            testCase.verifyFalse(contains(name,'.'));
        end

        function testEncodeNumber(testCase)
            x = [0 2 0.3 -1 1e-5 -1e-5 100000 0.49999 0.5];
            parts = arrayfun(@fem4mResult.encodeNumber,x,'UniformOutput',false);
            testCase.verifyEqual(parts,{'0','2','0p3','m1','1em05','m1em05', ...
                '100000','0p49999','0p5'});
            testCase.verifyEqual(cellfun(@fem4mResult.decodeNumber,parts),x);
            testCase.verifyFalse(any(contains(parts,{'.','-','_'})));
        end

        function testDecodeText(testCase)
            texts = {'p1','00','P2','top','1e5','-1','residual'};
            testCase.verifyTrue(all(isnan(cellfun(@fem4mResult.decodeNumber,texts))));
        end

        function testDisplayParts(testCase)
            parts = {'Lshape','Pk','uniform','residual','direct','2','0p3','150000'
                'Lshape','Pk','uniform','residual','direct','p1','m1em05','1000000'
                'Lshape','Pk','uniform','residual','direct','top','0','12345'
                'Lshape','Pk','uniform','residual','direct','m1','0','0'};
            values = fem4mResult.displayParts(parts);
            testCase.verifyEqual(values(:,1:5),parts(:,1:5));
            testCase.verifyEqual(values(:,6),{'2';'p1';'top';'-1'});
            testCase.verifyEqual(values(:,7),{'0.3';'-1e-05';'0';'0'});
            testCase.verifyEqual(values(:,8),{'15e4';'1e6';'12345';'0'});
        end

        function testDisplayLookup(testCase)
            % plotApp finds the file name of the shown values in its table
            names = {'Square_Pk_uniform_residual_direct_1_0p3_150000'
                'Square_Pk_uniform_residual_direct_1_0p49999_150000'
                'Square_Pk_uniform_residual_direct_1_0p5_150000'
                'Square_Pk_uniform_residual_direct_1_0p5_100000'};
            parts = cellfun(@(s)strsplit(s,'_'),names,'UniformOutput',false);
            values = fem4mResult.displayParts(vertcat(parts{:}));
            for k = 1 : numel(names)
                row = find(all(strcmp(values,repmat(values(k,:),numel(names),1)),2));
                testCase.verifyEqual(row,k);
            end
        end

        function testDefaultLabels(testCase)
            r = fem4mResult();
            testCase.verifyEqual(r.label,struct('problem','problem', ...
                'fe','fe','refine','refine','estimator','estimator', ...
                'solver','solver','param1','','param2','', ...
                'maxNrDoF','maxNrDoF'));
        end

        function testLabelMerge(testCase)
            r = fem4mResult(label=struct('problem','domain','param1','degree'));
            testCase.verifyEqual(r.label.problem,'domain');
            testCase.verifyEqual(r.label.param1,'degree');
            testCase.verifyEqual(r.label.fe,'fe');
            testCase.verifyEqual(r.label.param2,'');
        end

        function testSaveNameValue(testCase)
            levels = testCase.Levels;
            r = fem4mResultTest.pkResult();
            r.param2 = 0.5;
            r.label.param1 = 'degree';
            r.label.param2 = 'alpha';
            file = r.save(testCase.Folder,nrDoF=levels.nrDoF, ...
                H1Error=levels.H1Error,L2Error=levels.L2Error, ...
                ratio=levels.L2Error./levels.H1Error);
            testCase.verifyEqual(file,fullfile(testCase.Folder, ...
                'Square_Pk_uniform_residual_direct_2_0p5_10000.mat'));
            p = load(file).p;
            testCase.verifyEqual(p.problem,'Square');
            testCase.verifyEqual(p.fe,'Pk');
            testCase.verifyEqual(p.param1,2);
            testCase.verifyEqual(p.param2,0.5);
            testCase.verifyEqual(p.maxNrDoF,1e4);
            testCase.verifyEqual(p.label,{'problem','fe','refine', ...
                'estimator','solver','degree','alpha','maxNrDoF'});
            testCase.verifyNumElements(p.level,3);
            testCase.verifyEqual([p.level.nrDoF],levels.nrDoF);
            testCase.verifyEqual([p.level.H1Error],levels.H1Error);
            testCase.verifyEqual([p.level.L2Error],levels.L2Error);
            testCase.verifyEqual([p.level.ratio],levels.L2Error./levels.H1Error);
        end

        function testSaveStruct(testCase)
            r = fem4mResultTest.pkResult();
            levels = testCase.Levels;
            p1 = load(r.save(testCase.Folder,levels)).p;
            p2 = load(r.save(testCase.Folder,nrDoF=levels.nrDoF, ...
                H1Error=levels.H1Error,L2Error=levels.L2Error)).p;
            testCase.verifyEqual(p1,p2);
            testCase.verifyEqual(fieldnames(p1.level),{'nrDoF';'H1Error';'L2Error'});
        end

        function testSaveCreatesFolder(testCase)
            folder = fullfile(testCase.Folder,'Pk');
            file = fem4mResultTest.pkResult().save(folder,testCase.Levels);
            testCase.verifyTrue(isfile(file));
        end

        function testModuleResult(testCase)
            module = fem4mModuleElasticity();
            r = module.result(1e4,1,1,1,1,1);
            testCase.verifyEqual(r.param1,module.E);
            testCase.verifyEqual(r.param2,module.nu);
            testCase.verifyEqual({r.label.param1,r.label.param2},{'E','nu'});
            name = module.getFileName(1e4,1,1,1,1,1);
            testCase.verifyEqual(name,r.fileName());
            parts = strsplit(name,'_');
            testCase.verifyEqual(parts(6:8),{'100000','0p3','10000'});
        end

        function testModuleResultWithoutParameters(testCase)
            module = fem4mModuleElliptic();
            r = module.result(1e4,1,1,1,1,1);
            testCase.verifyEqual(r.fileName(),[module.problems{1},'_', ...
                module.finiteElements{1},'_uniform_residual_direct_0_0_10000']);
            testCase.verifyEqual(r.label.param1,'');
        end

    end

    methods (Static)

        function r = pkResult()
            %pkResult - Results file description of P2 with all slots set
            r = fem4mResult(problem='Square',fe='Pk',refine='uniform', ...
                estimator='residual',solver='direct',param1=2,param2=0, ...
                maxNrDoF=1e4);
        end

    end

end
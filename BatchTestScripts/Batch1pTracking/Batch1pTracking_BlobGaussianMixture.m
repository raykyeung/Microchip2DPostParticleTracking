%% Batch1pTracking33_pred_Final.m - [Script] Batch single p-Chip tracking
%
% Author: Raymond Yeung
% Release date: 2025
% E-mail: ryeun003@ucr.edu
% B2K Group, Dept. of Bioengineering, Univ. of California, Riverside
% Victor G. J. Rodgers Dept. of Bioengineering, Univ. of California, Riverside
% William H. Grover, Dept. of Bioengineering, Univ. of California, Riverside
% Philip L. Brisk Dept. of Computer Science, Univ. of California, Riverside

%% Start-up commands
close all force; %Close all figures
clear all; %Clear all variables, functions, scripts
% clearvars; %Clear variables
clc; %Clear the command window

profile clear; %Stop the Profiler and clear the recorded stats

%% Profiling code
%--To check the historysize needed, run/profile code, load profileHandle, run: Nrecorded = numel(profileHandle.FunctionHistory)
profile on -historysize 100000000

%% Timing - whole script, start
tScriptStart = tic; % pair 1: tic

%% Choose User
% User: UserOne
userOption.User = 'UserOne';

%% User Options (storing all variables that user may change to tweak code)
switch userOption.User
    case 'UserOne'
        % [SETTINGS: Video Summarization]
        % [Setting_01] Added frames before and after endpoint frames in which particle is detected
        userOption.addedFrames = 25;

        % [SETTINGS: Program modes]
        % [Setting_02] Drive priority: external, internal   (legacy drive search only)
        % [Setting_03] Processing modes = batch, single
        % [Setting_04] Run modes: debug, final
        % [Setting_05] Display modes: runWithDisplay, runWithoutDisplay
        % [Setting_06] Save location: defaultVideoLocation, anotherLocation
        % [Setting_07] Time function modes: timeFunction, noTimeFunction
        % [Setting_08] Time function data saving modes: timeSaveFunctionData, noTimeSaveFunctionData
        % [Setting_09] Parameter calculation: calcParameters, noCalcParameters

        userOption.drivePriority = 'internal';
        userOption.processingMode = 'batch';
        userOption.runMode = 'final';
        userOption.displayMode = 'runWithDisplay';
        userOption.saveLocation = 'defaultVideoLocation';
        userOption.evalTimeFunction = 'timeFunction';
        userOption.evalTimeSaveFunctionData = 'timeSaveFunctionData';
        userOption.evalCalcParameters = 'noCalcParameters';

        % Preferred data source: 'project' (recommended) or 'drives' (legacy)
        userOption.dataSourceMode = 'project';

        % [SETTINGS (OPTIONAL): Stored Data Path]
        % Keep empty for portable packaging
        userOption.dataPathOverride = '';

        % [SETTINGS (OPTIONAL): Drive]  (legacy search; keep empty for portability)
        userOption.driveOverride = '';

        % [SETTINGS (OPTIONAL): Single video analysis]
        userOption.singleVidOverride = '';

    otherwise
        error('Unknown user name.');
end

%% Initialize myFigures - Manages figure numbers and handles definition using structures containing fields with dynamic expression
myFigures = struct;

myFigures.fig = struct;
myFigures.handle = strings(1,1);
myFigures.figNum = 0;
myFigures.figNumArray = [];
myFigures.counterNum = 0;
myFigures.counterNumArray = [];
myFigures.processNum = 0;
myFigures.processNumArray = [];
myFigures.initFigNumLoop = 0;
myFigures.flagInLoop = false;

%% Flags
[settingAddedFrames,flagInternalDrive,flagBatch,flagDebug,flagRunWithDisplay,flagSaveLocDefault,flagTimeFunction,flagTimeSaveFunctionData,flagCalcParameters,flagDataPathOverride,flagDriveOverride,flagSingleVidOverride] = ...
    B2KFlagsBGM( ...
        userOption.addedFrames, ...
        userOption.drivePriority, ...
        userOption.processingMode, ...
        userOption.runMode, ...
        userOption.displayMode, ...
        userOption.saveLocation, ...
        userOption.evalTimeFunction, ...
        userOption.evalTimeSaveFunctionData, ...
        userOption.evalCalcParameters, ...
        userOption.dataPathOverride, ...
        userOption.driveOverride, ...
        userOption.singleVidOverride);

% Local flag for project-root data source (for packaged project)
flagProjectSource = strcmpi(userOption.dataSourceMode, 'project');

%% Find folder storing data videos and set up output folder

% Saved default drives / data videos folder / output folder locations
ExpInternalDriveNames = ["INT_DRIVE_A","INT_DRIVE_B","INT_DRIVE_C"];
ExpExternalDriveNames = ["EXT_DRIVE_A","EXT_DRIVE_B","EXT_DRIVE_C"];

% List all physical drives on computer (used only in legacy mode)
infoTable = B2KlistPhysicalDrives;

% Warn if both overrides are present
if flagDataPathOverride && flagDriveOverride
    warning('Both dataPathOverride and driveOverride provided. Defaulting to dataPathOverride ...');
end

if flagDataPathOverride
    % 1) Hard override of data path takes precedence
    selpath = userOption.dataPathOverride;
    fprintf('Analyzing data from chosen (overridden) location, %s ...\n', selpath);

elseif flagProjectSource
    % 2) Project-root anchored (recommended for packaged projects)
    projRoot = getProjectRootOrHere();  % absolute path to the MATLAB Project root (or this file’s folder)
    dv = fullfile(projRoot, 'Data_Videos');  % videos expected in <project root>/Data_Videos/*.avi
    if exist(dv, 'dir') ~= 7
        warning('[B2K] Expected "<project>/Data_Videos" not found. Please select the working folder that contains "Data_Videos".');
        selpath = uigetdir;
        if isequal(selpath,0), error('No folder selected.'); end
        fprintf('Analyzing data from chosen (GUI) location, %s ...\n', selpath);
    else
        selpath = projRoot; % anchor to project root for reading videos
        fprintf('Analyzing data from Project root: %s\n', selpath);
    end

else
    % 3) Legacy drive search (compatible, but not portable)
    if ~flagDriveOverride
        if flagInternalDrive
            ExpDriveNames = [ExpInternalDriveNames, ExpExternalDriveNames];
        else
            ExpDriveNames = [ExpExternalDriveNames, ExpInternalDriveNames];
        end

        matchDriveName = [];
        for d = 1:numel(ExpDriveNames)
            if ispc
                matchDriveName = infoTable.DeviceID(strcmp(infoTable.VolumeName, ExpDriveNames(d)), :);
            elseif ismac
                matchDriveName = fullfile('/Volumes', infoTable.VolumeName(strcmp(infoTable.VolumeName, ExpDriveNames(d)), :));
            else % if isunix
                error('Have not tested in Linux');
            end
            if ~isempty(matchDriveName), break; end
        end

        if ~isempty(matchDriveName)
            selpath = matchDriveName;
            if flagInternalDrive
                if d <= numel(ExpInternalDriveNames)
                    fprintf('Analyzing data from internal drive, %s%s ...\n', matchDriveName, ExpDriveNames(d));
                else
                    fprintf('Analyzing data from external drive, %s%s ...\n', matchDriveName, ExpDriveNames(d));
                end
            else
                if d > numel(ExpInternalDriveNames)
                    fprintf('Analyzing data from internal drive, %s%s ...\n', matchDriveName, ExpDriveNames(d));
                else
                    fprintf('Analyzing data from external drive, %s%s ...\n', matchDriveName, ExpDriveNames(d));
                end
            end
        else
            warning('[B2K] Default drives not found. Please select the working folder. The working folder should have a subfolder titled "Data_Videos" containing the experimental videos.');
            selpath = uigetdir; % absolute path to selected working folder
            if isequal(selpath,0), error('No folder selected.'); end
            fprintf('Analyzing data from chosen (GUI) location, %s ...\n', selpath);
        end
    else
        selpath = userOption.driveOverride;
        fprintf('Analyzing data from chosen (overridden) drive, %s ...\n', selpath);
    end
end

%% Setup main output folder and output subfolders
% For project-root mode, write outputs to OS "Downloads/OutputFolder"; else write under selpath.
outputBase = selpath;
if flagProjectSource
    outputBase = getDownloadsFolder();   % OS-appropriate Downloads folder
    fprintf('Output base set to Downloads folder: %s\n', outputBase);
end

%[Folder structure expectation for packaged projects]
% <project root>
% ├─ Data_Videos
% │    └─ *.avi
% └─ (Outputs go to) <Downloads>/OutputFolder/<Child>/<Child_yyMMdd_xx>/...

%5a Current script name ('Batch1pTracking01') with robust fallback
currScriptName = currentScriptNameOr('MainScript');

%5b Date ('yyMMdd')
currDate = char(datetime('now','Format','yyMMdd'));

%5 Combine current script name and date ('Batch1pTracking01_230613')
currScriptNameAddDate = append(currScriptName,'_',currDate);

%4 Child output folder ('Batch1pTracking') extracted from script name
childOutputFolderName = char(regexp(currScriptName,'[a-zA-Z]+\w*[a-zA-Z]+','match','once')); % first match only
if isempty(childOutputFolderName), childOutputFolderName = currScriptName; end

%3 Parent output folder name (container under outputBase)
parentOutputFolderName = 'OutputFolder';

%4 Absolute file path for child output folder
childOutputFolderPath = fullfile(outputBase, parentOutputFolderName, childOutputFolderName);

%5 Absolute file path for grandchild output folder
grandchildOutputFolderUnnumberedPath = fullfile(childOutputFolderPath, currScriptNameAddDate);
grandchildOutputFolderNumberedName   = B2Knextname(grandchildOutputFolderUnnumberedPath,'_01','');
grandchildOutputFolderPath           = fullfile(outputBase, parentOutputFolderName, childOutputFolderName, grandchildOutputFolderNumberedName);

% Ensure folders exist
mkIfMissing(childOutputFolderPath);
mkIfMissing(grandchildOutputFolderPath);

% Video summarization output folder path
vidSumOutputFolderPath = fullfile(grandchildOutputFolderPath,'VideoSummarization_Output');
mkIfMissing(vidSumOutputFolderPath);

% Channel wall detection output folder path
wallDetOutputFolderPath = fullfile(grandchildOutputFolderPath,'ChannelWallDetection_Output');
mkIfMissing(wallDetOutputFolderPath);

% 1p Tracking output folder path
onePTrackOutputFolderPath = fullfile(grandchildOutputFolderPath,'OnePTracking_Output');
mkIfMissing(onePTrackOutputFolderPath);

%% Create diary file for logging commands, keyboard inputs, text output
diaryOutputName       = append(grandchildOutputFolderNumberedName,'_diary');
diaryOutputFolderPath = fullfile(grandchildOutputFolderPath, diaryOutputName);
diary(diaryOutputFolderPath)

%% GUI functions
% uigetfile, uiputfile, uigetdir, uiopen, uisave

%% Video datastore
% Project mode = use ONLY "<project root>/Data_Videos"
% Legacy mode   = use "Data_Videos" (and optional _2/_3) under selpath
if flagProjectSource
    baseVideos = fullfile(selpath, 'Data_Videos');  % <project root>/Data_Videos
    if exist(baseVideos,'dir') ~= 7
        error('Project mode: folder "%s" not found. Expected at "<project root>/Data_Videos".', baseVideos);
    end

    if flagBatch
        % Batch: recurse inside Data_Videos only
        fds = fileDatastore(baseVideos, ...
            'ReadFcn', @myread, ...
            'IncludeSubfolders', true, ...
            'FileExtensions', '.avi');
        vidPathStack = fds.Files;
        totalNumVideosToAnalyze = numel(vidPathStack);
        if totalNumVideosToAnalyze == 0
            error('No .avi videos found under "%s".', baseVideos);
        end
    else
        % Single: pick a file starting from Data_Videos
        totalNumVideosToAnalyze = 1;
        if flagSingleVidOverride && ~isempty(userOption.singleVidOverride)
            warning('Single video mode: using provided singleVidOverride ...');
            vidPathStack = userOption.singleVidOverride;
        else
            warning('Single video mode: please select a video from the "Data_Videos" folder.');
            figSingleVidAlert = uifigure('WindowStyle','alwaysontop');

            monitorInfo = get(0,"MonitorPositions");
            [monitorResScaled,~] = max(monitorInfo(:,3));
            maxIdx = find(monitorInfo==monitorResScaled);
            szMonitorInfo = size(monitorInfo);
            [monitorNum,~] = ind2sub(szMonitorInfo,maxIdx);
            if size(monitorNum,1) > 1
                monitorChoice = min(monitorNum);
            else
                monitorChoice = monitorNum;
            end
            figSingleVidAlert.Position = [monitorInfo(monitorChoice,1),monitorInfo(monitorChoice,2),450,175];
            movegui(figSingleVidAlert,"center");

            uialert(figSingleVidAlert, ...
                'Program currently set to analyze a single video. Please select a video from the "Data_Videos" folder.', ...
                'Single video processing mode', ...
                'Modal', false, 'Icon', 'info', 'CloseFcn', 'uiresume(figSingleVidAlert)');
            uiwait(figSingleVidAlert);
            close(figSingleVidAlert);

            [vidName,vidPath] = uigetfile({'*.avi','AVI Video (*.avi)'}, 'Select video', baseVideos);
            if isequal(vidName,0)
                error('No video selected.');
            end
            vidPathStack = fullfile(vidPath,vidName);
        end
    end

else
    % LEGACY mode (kept compatible). Allow Data_Videos_2/_3 but only if they exist.
    loc1 = fullfile(selpath,'Data_Videos');
    loc2 = fullfile(selpath,'Data_Videos_2');
    loc3 = fullfile(selpath,'Data_Videos_3');
    locAll = {loc1,loc2,loc3};
    locAll = locAll(cellfun(@(p) exist(p,'dir')==7, locAll));  % keep existing only
    locAll = ensureCharCell(locAll);  % make sure it's a cell array of char (for older MATLAB)

    if flagBatch
        if isempty(locAll)
            error('No "Data_Videos" folders found under %s. Expected at least "<root>/Data_Videos".', selpath);
        end
        fds = fileDatastore(locAll, ...
            'ReadFcn', @myread, ...
            'IncludeSubfolders', true, ...
            'FileExtensions', '.avi');
        vidPathStack = fds.Files;
        totalNumVideosToAnalyze = numel(vidPathStack);
        if totalNumVideosToAnalyze == 0
            error('No .avi videos found under %s.', strjoin(locAll, '; '));
        end
    else
        totalNumVideosToAnalyze = 1;
        if flagSingleVidOverride && ~isempty(userOption.singleVidOverride)
            warning('Single video mode: using provided singleVidOverride ...');
            vidPathStack = userOption.singleVidOverride;
        else
            warning('Single video mode: please select a video from the "Data_Videos" folder.');
            figSingleVidAlert = uifigure('WindowStyle','alwaysontop');

            monitorInfo = get(0,"MonitorPositions");
            [monitorResScaled,~] = max(monitorInfo(:,3));
            maxIdx = find(monitorInfo==monitorResScaled);
            szMonitorInfo = size(monitorInfo);
            [monitorNum,~] = ind2sub(szMonitorInfo,maxIdx);
            if size(monitorNum,1) > 1
                monitorChoice = min(monitorNum);
            else
                monitorChoice = monitorNum;
            end
            figSingleVidAlert.Position = [monitorInfo(monitorChoice,1),monitorInfo(monitorChoice,2),450,175];
            movegui(figSingleVidAlert,"center");

            uialert(figSingleVidAlert, ...
                'Program currently set to analyze a single video. Please select a video from the "Data_Videos" folder.', ...
                'Single video processing mode', ...
                'Modal', false, 'Icon', 'info', 'CloseFcn', 'uiresume(figSingleVidAlert)');
            uiwait(figSingleVidAlert);
            close(figSingleVidAlert);

            % Prefer starting in "<selpath>/Data_Videos" if it exists
            startDir = loc1; if exist(startDir,'dir') ~= 7, startDir = selpath; end
            [vidName,vidPath] = uigetfile({'*.avi','AVI Video (*.avi)'}, 'Select video', startDir);
            if isequal(vidName,0)
                error('No video selected.');
            end
            vidPathStack = fullfile(vidPath,vidName);
        end
    end
end

%% Extract filtered, ordered subset of data videos and create a new file datastore

if ~flagProjectSource
    fdsFilePaths = fds.Files;
    
    [~, names, exts] = cellfun(@fileparts, fdsFilePaths, 'UniformOutput', false);
    fullNames = strcat(names, exts);
    
    orderedSubset = filteredTableF.FileName;
    
    orderedPaths = cell(size(orderedSubset));
    for k = 1:numel(orderedSubset)
      idx = find(strcmp(fullNames, orderedSubset{k}), 1);
      if isempty(idx)
        error('File "%s" not found in original datastore.', orderedSubset{k});
      end
      orderedPaths{k} = fdsFilePaths{idx};
    end
    
    filteredFds = fileDatastore(orderedPaths,'ReadFcn',@myread,'FileExtensions','.avi');
else
    filteredFds = fds;
end

%% Video Summarization setup
addedFrames = userOption.addedFrames;
idx_BatchVideoSummarization = 14;
cell_BatchVideoSummarization = cell(size(10,1),idx_BatchVideoSummarization);
vidSumPassCount = 0;
vidSumErrCount = 0;

%% Channel Wall Detection setup
idx_cell_BatchChannelDetection = 15;
cell_BatchChannelDetection = cell(size(10,1),idx_cell_BatchChannelDetection);

cell_AllData = cell(2,2);

HLcorrectionidx = 0;

%% Particle Tracking setup
idx_cell_BatchParticleTracking = 14;
cell_BatchParticleTracking = cell(size(10,1),idx_cell_BatchParticleTracking);

%% Loop

loopStart = 1;
loopEnd = 1;
% loopEnd = length(filteredFds.Files);

%% Timing function process and function data saving using tic and toc

if flagTimeFunction == true
    tEndFnVidSum = NaN(size(vidPathStack,1),1);
    tEndSaveFnVidSumData = NaN(size(vidPathStack,1),1);
    
    tEndFnCoordSys = NaN(size(vidPathStack,1),1);
    tEndSaveFnCoordSysData = NaN(size(vidPathStack,1),1);
    
    tEndFnParticleTrack = NaN(size(vidPathStack,1),1);
    tEndSaveFnParticleTrackData = NaN(size(vidPathStack,1),1);
end

for vidNumIdx = loopStart:loopEnd
    fprintf('Processing video #%d...video %d/%d...%s...\n', vidNumIdx,vidNumIdx-loopStart+1,loopEnd-loopStart+1,datetime)
    
    %% Update flag loop status
    myFigures.flagInLoop = true;

    if vidNumIdx == loopStart
        myFigures.initFigNumLoop = myFigures.figNum;
    end

    %% Read video to create video object and vid name
    if flagBatch == true
        vidPath = filteredFds.Files{vidNumIdx};
        vidObj = myread(vidPath);
        [vidFolder,vidNameNoExt,ext] = fileparts(vidPath);
        vid = append(vidNameNoExt,ext);
    else %if flagBatch == false
        vidPath = vidPathStack;
        vidObj = VideoReader(vidPath);
        [vidFolder,vidNameNoExt,ext] = fileparts(vidPath);
        vid = append(vidNameNoExt,ext);
    end

    %% B2KVideoSummary - Video summarization to identify start and end frames to analyze

    %% Timing start, FnVidSum
    if flagTimeFunction == true
        tStartFnVidSum = tic;
    end

    try
        [L_start,R_end,tFrames,nFrames,...
            flag_pChip_present,flag_pChip_stuck,flag_no_pChip,flag_video_corrupt,...
            myFigures] = B2KVideoSummaryBGM(vidObj,myFigures,flagDebug,flagRunWithDisplay,addedFrames);

        %% Timing end, FnVidSum
        if flagTimeFunction == true
            tEndFnVidSum(vidNumIdx,1) = toc(tStartFnVidSum);
        end

        %% Timing start, SaveFnVidSumData
        if flagTimeFunction == true
            tStartSaveFnVidSumData = tic;
        end
    catch ME
        %% Timing end, FnVidSum
        if flagTimeFunction == true
            tEndFnVidSum(vidNumIdx,1) = toc(tStartFnVidSum);
        end

        %% Timing start, SaveFnVidSumData
        if flagTimeFunction == true
            tStartSaveFnVidSumData = tic;
        end

        %% Manually add appropriate error message(s) from function call(s)
        % Construct a detailed error message including the stack trace
        errorMsg = sprintf('\n%s\n\n', ME.message);
    
        for k = 1:length(ME.stack)-1
            % Read the file and get the specific line where the error occurred
            fileContent = fileread(ME.stack(k).file);
            fileLines = strsplit(fileContent, '\n');
            errorLine = fileLines{ME.stack(k).line};
    
            % Create clickable bold link using MATLAB's hyperlink format with HTML bold tags
            errorMsg = sprintf('%sError in <a href="matlab: opentoline(''%s'',%d)">%s</a> (line %d):\n%s\n', ...
                               errorMsg, ME.stack(k).file, ME.stack(k).line, ME.stack(k).name, ME.stack(k).line, errorLine);
        end
    
        % Use the error function to display the error message in red and stop execution
        error(errorMsg);

        %%
        % -- If error in B2KVideoSummary, then: --
        % 1) Save table results to .mat
        % 2) Save used functions to log .txt

        table_BatchVideoSummarization = cell2table(cell_BatchVideoSummarization,...
            "VariableNames",...
            ["Video Reader Object" ...                      % 1
            "Video File Name" ...                           % 2
            "Time Video Summary (s)" ...                    % 3
            "Time Save Video Summary (s)" ...               % 4
            "Error Count" ...                               % 5
            "Pass Count" ...                                % 6
            "Start Frame" ...                               % 7    
            "End Frame" ...                                 % 8
            "Number of Frames to Analyze" ...               % 9
            "Total Number of Frames" ...                    % 10
            "Flag p-Chip Present" ...                       % 11
            "Flag p-Chip Stuck" ...                         % 12
            "Flag no p-Chip" ...                            % 13
            "Flag Video Corrupt"]);                         % 14

        %% 1) Save table results to .mat
        vidSumTableOutputName = append("table_",grandchildOutputFolderNumberedName,"_vidSum",".mat");
        vidSumTableOutputPath = fullfile(vidSumOutputFolderPath,vidSumTableOutputName);

        save(vidSumTableOutputPath,"table_BatchVideoSummarization");

        %% 2) Save used functions to log .txt file
        currScriptNameWithExt = append(currScriptName,'.m');
        [fList,pList] = matlab.codetools.requiredFilesAndProducts(currScriptNameWithExt);

        logFileName = append(grandchildOutputFolderNumberedName,'_log','.txt');
        logFileFullPath = fullfile(grandchildOutputFolderPath,logFileName);
        writelines([newline,'Function run (function_YYMMDD_instance):'],logFileFullPath);
        writelines(grandchildOutputFolderNumberedName,logFileFullPath,"WriteMode","append");
        
        writelines([newline,'Main function:'],logFileFullPath,"WriteMode","append");
        writelines(currScriptNameWithExt,logFileFullPath,"WriteMode","append");

        writelines([newline,'Settings:'],logFileFullPath,"WriteMode","append");
        writelines(sprintf('Added frames: %d',userOption.addedFrames),logFileFullPath,"WriteMode","append");
        writelines(sprintf('Video number start: %d',loopStart),logFileFullPath,"WriteMode","append");
        writelines(sprintf('Video number end: %d',loopEnd),logFileFullPath,"WriteMode","append");

        writelines([newline,'Flags: '],logFileFullPath,"WriteMode","append");
        writelines(['Drive priority: ',userOption.drivePriority],logFileFullPath,"WriteMode","append");
        writelines(['Processing mode: ',userOption.processingMode],logFileFullPath,"WriteMode","append");
        writelines(['Run mode: ',userOption.runMode],logFileFullPath,"WriteMode","append");
        writelines(['Display mode: ',userOption.displayMode],logFileFullPath,"WriteMode","append");
        writelines(['Save location: ',userOption.saveLocation],logFileFullPath,"WriteMode","append");
        writelines(['Time process function: ',userOption.evalTimeFunction],logFileFullPath,"WriteMode","append");
        writelines(['Time saving process function data: ',userOption.evalTimeSaveFunctionData],logFileFullPath,"WriteMode","append");
        writelines(['Data path override: ',userOption.dataPathOverride],logFileFullPath,"WriteMode","append");
        writelines(['Drive override: ',userOption.driveOverride],logFileFullPath,"WriteMode","append");
        writelines(['Single video override: ',userOption.singleVidOverride],logFileFullPath,"WriteMode","append");

        writelines([newline,'Functions used:'],logFileFullPath,"WriteMode","append");
        for idxfList = 1:length(fList)
            fListFullPath = fList{idxfList};
            [fListFilePath,fListName,FListExt] = fileparts(fListFullPath);
            writelines(fListName,logFileFullPath,"WriteMode","append");
        end

        writelines([newline,'Required program files:'],logFileFullPath,"WriteMode","append");

        for idxpList = 1:length(pList)
            pListName = pList(idxpList).Name;
            pListVer = pList(idxpList).Version;
            pListNameVer = append(pListName,"_",pListVer);
            writelines(pListNameVer,logFileFullPath,"WriteMode","append")
        end

        %% Disable diary logging
        diary off

        %% Clear figure every loop iteration; close at end of loop
        if vidNumIdx == loopEnd
            for idxFig = 1:myFigures.processNum
                idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
                close(myFigures.fig.(myFigures.handle(idxFigRef)))
            end
            myFigures.flagInLoop = false;
        else
            for idxFig = 1:myFigures.processNum
                idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
                clf(myFigures.fig.(myFigures.handle(idxFigRef)))
            end
        end

        %% Throw default error message(s)
        rethrow(ME)
    end

    % -- If no error in B2KVideoSummary, then: --
    % 3) Save workspace variables in .mat
    % 4) Save figure in .tif
    % 5) Store data in cell array

    %% 3) Save workspace variables for each individual run in .mat file
    vidSumFileName = append(grandchildOutputFolderNumberedName,'_vidSum');
    vidSumFileNameUnnumberedPath = fullfile(vidSumOutputFolderPath,vidSumFileName);
    vidSumFileName = B2Knextname(vidSumFileNameUnnumberedPath,'_vid_01','.mat');

    fileNamePath = fullfile(vidSumOutputFolderPath,vidSumFileName);

    % Exclude any figures or graphics objects while allowing to keep open
    % Get a list of all variables
    allvars = whos;

    % Identify the variables that ARE NOT graphics handles. This uses a regular
    % expression on the class of each variable to check if it's a graphics object
    tosave1 = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));

    % Manually remove 'myFigures' struct containing figure handles
    tosave2 = cellfun(@isempty, regexp({allvars.name},'myFigures'));

    tosavefinal = tosave1 & tosave2;

    % Pass these variable names to save
    save(fileNamePath, allvars(tosavefinal).name,'-v7.3');

    %% 4) Save figure

    for idxFig = 1:myFigures.processNum
        idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig; %[250312 - ERROR, -2 value for second vid analyzed]

        % Export figure as .tif (for publication)
        fileImg = strrep(fileNamePath,'.mat','.tif');
        figNumStr = sprintf('_fig%d',idxFigRef);
        fileNameImg = insertBefore(fileImg,'.tif',figNumStr);
        % exportgraphics(myFigures.fig.(myFigures.handle(myFigures.figNum)),fileNameImg1,'Resolution',600);
        print(myFigures.fig.(myFigures.handle(idxFigRef)),fileNameImg,'-dtiffn');

        % Set CreateFcn callback
        set(myFigures.fig.(myFigures.handle(myFigures.figNum)), 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 

        % Save figure as .mat file (for reference and debugging)
        fileFig = strrep(fileNamePath,'.mat','.fig');
        fileNameFig = insertBefore(fileFig,'.fig',figNumStr);
        savefig(myFigures.fig.(myFigures.handle(idxFigRef)),fileNameFig);
    end

    %% 5) Store data in cell array
    cell_BatchVideoSummarization{vidNumIdx,1} = vidObj;
    cell_BatchVideoSummarization{vidNumIdx,2} = vidObj.Name;
    cell_BatchVideoSummarization{vidNumIdx,3} = tEndFnVidSum(vidNumIdx,1);
    cell_BatchVideoSummarization{vidNumIdx,4} = NaN;
    cell_BatchVideoSummarization{vidNumIdx,5} = NaN;
    cell_BatchVideoSummarization{vidNumIdx,6} = NaN;
    cell_BatchVideoSummarization{vidNumIdx,7} = L_start;
    cell_BatchVideoSummarization{vidNumIdx,8} = R_end;
    cell_BatchVideoSummarization{vidNumIdx,9} = tFrames;
    cell_BatchVideoSummarization{vidNumIdx,10} = nFrames;
    cell_BatchVideoSummarization{vidNumIdx,11} = flag_pChip_present;
    cell_BatchVideoSummarization{vidNumIdx,12} = flag_pChip_stuck;
    cell_BatchVideoSummarization{vidNumIdx,13} = flag_no_pChip;
    cell_BatchVideoSummarization{vidNumIdx,14} = flag_video_corrupt;

    %%
    if isnan(flag_pChip_present) || ~isnan(flag_pChip_stuck) || ~isnan(flag_no_pChip) || ~isnan(flag_video_corrupt)
        vidSumErrCount = vidSumErrCount + 1;
        cell_BatchVideoSummarization{vidNumIdx,5} = vidSumErrCount;
        clf(myFigures.fig.(myFigures.handle(myFigures.figNum)))

        %% 3) Save workspace variables for each individual run in .mat file
        wallDetFileName = append(grandchildOutputFolderNumberedName,'_wallDet');
        wallDetFileNameUnnumberedPath = fullfile(wallDetOutputFolderPath,wallDetFileName);
        wallDetFileName = B2Knextname(wallDetFileNameUnnumberedPath,'_vid_01','.mat');

        fileNamePath = fullfile(wallDetOutputFolderPath,wallDetFileName);

        % Exclude any figures or graphics objects while allowing to keep open
        % Get a list of all variables
        allvars = whos;

        % Identify the variables that ARE NOT graphics handles. This uses a regular
        % expression on the class of each variable to check if it's a graphics object
        tosave1 = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));

        % Manually remove 'myFigures' struct containing figure handles
        tosave2 = cellfun(@isempty, regexp({allvars.name},'myFigures'));

        tosavefinal = tosave1 & tosave2;

        % Pass these variable names to save
        save(fileNamePath, allvars(tosavefinal).name,'-v7.3');

        %% 4) Save figure

        for idxFig = 1:myFigures.processNum
            idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;

            % Export figure as .tif (for publication)
            fileImg = strrep(fileNamePath,'.mat','.tif');
            figNumStr = sprintf('_fig%d',idxFigRef);
            fileNameImg = insertBefore(fileImg,'.tif',figNumStr);
            % exportgraphics(myFigures.fig.(myFigures.handle(myFigures.figNum)),fileNameImg1,'Resolution',600);
            print(myFigures.fig.(myFigures.handle(idxFigRef)),fileNameImg,'-dtiffn');

            % Set CreateFcn callback
            set(myFigures.fig.(myFigures.handle(myFigures.figNum)), 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 

            % Save figure as .mat file (for reference and debugging)
            fileFig = strrep(fileNamePath,'.mat','.fig');
            fileNameFig = insertBefore(fileFig,'.fig',figNumStr);
            savefig(myFigures.fig.(myFigures.handle(idxFigRef)),fileNameFig);
        end

        %% 5) Store data in cell array
        cell_BatchChannelDetection{vidNumIdx,1} = myread(filteredFds.Files{vidNumIdx});
        cell_BatchChannelDetection{vidNumIdx,2} = cell_BatchChannelDetection{vidNumIdx,1}.Name;

        for j = 3:idx_cell_BatchChannelDetection
            cell_BatchChannelDetection{vidNumIdx,j} = NaN;
        end

        %% Clear figure every loop iteration; close at end of loop
        if vidNumIdx == loopEnd
            for idxFig = 1:myFigures.processNum
                idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
                close(myFigures.fig.(myFigures.handle(idxFigRef)))
            end
            myFigures.flagInLoop = false;
        else
            for idxFig = 1:myFigures.processNum
                idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
                clf(myFigures.fig.(myFigures.handle(idxFigRef)))
            end
        end

        continue % pass control to next iteration of for loop since an experimental error has been flagged

    end

    %%
    vidSumPassCount = vidSumPassCount + 1;
    cell_BatchVideoSummarization{vidNumIdx,6} = vidSumPassCount;

    %% Clear figure every loop iteration; close at end of loop
    if vidNumIdx == loopEnd
        for idxFig = 1:myFigures.processNum
            idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
            close(myFigures.fig.(myFigures.handle(idxFigRef)))
        end
        myFigures.flagInLoop = false;
    else
        for idxFig = 1:myFigures.processNum
            idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
            clf(myFigures.fig.(myFigures.handle(idxFigRef)))
        end
    end

    %% Timing end, SaveFnVidSumData
    if flagTimeFunction == true
        tEndSaveFnVidSumData(vidNumIdx,1) = toc(tStartSaveFnVidSumData);
    end

    cell_BatchVideoSummarization{vidNumIdx,4} = tEndSaveFnVidSumData(vidNumIdx,1);

    %% Read results from BatchVideoSummarization and skip videos with errors
    if isnan(cell_BatchVideoSummarization{vidNumIdx,6})

        cell_BatchChannelDetection{vidNumIdx,1} = myread(filteredFds.Files{vidNumIdx});
        cell_BatchChannelDetection{vidNumIdx,2} = cell_BatchChannelDetection{vidNumIdx,1}.Name;

        for j = 3:idx_cell_BatchChannelDetection
            cell_BatchChannelDetection{vidNumIdx,j} = NaN;
        end

        continue
    end

    %% Read video summarization data
    L_start = cell_BatchVideoSummarization{vidNumIdx,7};
    R_end = cell_BatchVideoSummarization{vidNumIdx,8};
    tFrames = cell_BatchVideoSummarization{vidNumIdx,9};

    L_start_orig = L_start;
    R_end_orig = R_end;
    tFrames_orig = tFrames;

    %% Experimental variables
    [EXP_numpChips,EXP_solvent,EXP_run,EXP_date,EXP_H_channel,EXP_W_channel,EXP_flowrate,EXP_desiredframerate,EXP_collection] = B2KExpParam(vid);
    EXP_actualframerate = vidObj.FrameRate;

    %% B2KCoordinateSystem - Image calibration involving channel wall detection and mapping data to appropriate coordinate system
    CS_string = 'CR';

    %% Timing start, FnCoordSys
    if flagTimeFunction == true
        tStartFnCoordSys = tic;
    end

    try
        [I, CS, WR, CR, NR, x_pos_caption, y_pos_caption,myFigures,...
            resizeScale, edgeMethod, R_highest_resolution,T_highest_resolution,...
            H_highest_numpeaks,Hlines_highest_fillgap,Hlines_highest_minlength,...
            Hlines_highest,Hlines_highest_num] = B2KCoordinateSystemBGM(vidObj,L_start_orig, HLcorrectionidx,EXP_W_channel,CS_string,myFigures,flagDebug,flagRunWithDisplay);
    
        %% Timing end, FnCoordSys
        if flagTimeFunction == true
            tEndFnCoordSys(vidNumIdx,1) = toc(tStartFnCoordSys);
        end

        %% Timing start, SaveFnCoordSysData
        if flagTimeFunction == true
            tStartSaveFnCoordSysData = tic;
        end

    catch ME
        %% Timing end, FnCoordSys
        if flagTimeFunction == true
            tEndFnCoordSys(vidNumIdx,1) = toc(tStartFnCoordSys);
        end

        %% Timing start, SaveFnCoordSysData
        if flagTimeFunction == true
            tStartSaveFnCoordSysData = tic;
        end

        %% Manually add appropriate error message(s) from function call(s)
        % Construct a detailed error message including the stack trace
        errorMsg = sprintf('\n%s\n\n', ME.message);
    
        for k = 1:length(ME.stack)-1
            % Read the file and get the specific line where the error occurred
            fileContent = fileread(ME.stack(k).file);
            fileLines = strsplit(fileContent, '\n');
            errorLine = fileLines{ME.stack(k).line};
    
            % Create clickable bold link using MATLAB's hyperlink format with HTML bold tags
            errorMsg = sprintf('%sError in <a href="matlab: opentoline(''%s'',%d)">%s</a> (line %d):\n%s\n', ...
                               errorMsg, ME.stack(k).file, ME.stack(k).line, ME.stack(k).name, ME.stack(k).line, errorLine);
        end
    
        % Use the error function to display the error message in red and stop execution
        error(errorMsg);

        %%
        % -- If error in B2KCoordinateSystemBatch, then: --
        % 1) Save table results to .mat
        % 2) Save used functions to log .txt

        table_BatchChannelDetection = cell2table(cell_BatchChannelDetection,...
            "VariableNames",...
            ["Video Reader Object" ...                      % 1
            "Video File Name" ...                           % 2
            "Time Coordinate System (s)"                    % 3
            "Time Save Coordinate System (s)"               % 4
            "Start Frame" ...                               % 5
            "End Frame" ...                                 % 6
            "Resize Scale" ...                              % 7
            "Edge Method" ...                               % 8
            "Hough Rho Resolution" ...                      % 9
            "Hough Theta Resolution" ...                    % 10
            "Hough Number of Peaks" ...                     % 11
            "Houghlines Fill Gap" ...                       % 12
            "Houghlines Minimum Length" ...                 % 13
            "Houghlines" ...                                % 14
            "Houghlines Number Detected"]);                 % 15

        %% 1) Save table results to .mat
        wallDetTableOutputName = append("table_",grandchildOutputFolderNumberedName,"_wallDet",".mat");
        wallDetTableOutputPath = fullfile(wallDetOutputFolderPath,wallDetTableOutputName);

        % save(fileNamePath);
        save(wallDetTableOutputPath,"table_BatchChannelDetection",'-v7.3');

        %% 2) Save used functions to log .txt file

        currScriptNameWithExt = append(currScriptName,'.m');
        [fList,pList] = matlab.codetools.requiredFilesAndProducts(currScriptNameWithExt);

        logFileName = append(grandchildOutputFolderNumberedName,'_log','.txt');
        logFileFullPath = fullfile(grandchildOutputFolderPath,logFileName);
        writelines([newline,'Function run (function_YYMMDD_instance):'],logFileFullPath);
        writelines(grandchildOutputFolderNumberedName,logFileFullPath,"WriteMode","append");
        
        writelines([newline,'Main function:'],logFileFullPath,"WriteMode","append");
        writelines(currScriptNameWithExt,logFileFullPath,"WriteMode","append");

        writelines([newline,'Settings:'],logFileFullPath,"WriteMode","append");
        writelines(sprintf('Added frames: %d',userOption.addedFrames),logFileFullPath,"WriteMode","append");
        writelines(sprintf('Video number start: %d',loopStart),logFileFullPath,"WriteMode","append");
        writelines(sprintf('Video number end: %d',loopEnd),logFileFullPath,"WriteMode","append");

        writelines([newline,'Flags: '],logFileFullPath,"WriteMode","append");
        writelines(['Drive priority: ',userOption.drivePriority],logFileFullPath,"WriteMode","append");
        writelines(['Processing mode: ',userOption.processingMode],logFileFullPath,"WriteMode","append");
        writelines(['Run mode: ',userOption.runMode],logFileFullPath,"WriteMode","append");
        writelines(['Display mode: ',userOption.displayMode],logFileFullPath,"WriteMode","append");
        writelines(['Save location: ',userOption.saveLocation],logFileFullPath,"WriteMode","append");
        writelines(['Time process function: ',userOption.evalTimeFunction],logFileFullPath,"WriteMode","append");
        writelines(['Time saving process function data: ',userOption.evalTimeSaveFunctionData],logFileFullPath,"WriteMode","append");
        writelines(['Data path override: ',userOption.dataPathOverride],logFileFullPath,"WriteMode","append");
        writelines(['Drive override: ',userOption.driveOverride],logFileFullPath,"WriteMode","append");
        writelines(['Single video override: ',userOption.singleVidOverride],logFileFullPath,"WriteMode","append");

        writelines([newline,'Functions used:'],logFileFullPath,"WriteMode","append");
        for idxfList = 1:length(fList)
            fListFullPath = fList{idxfList};
            [fListFilePath,fListName,FListExt] = fileparts(fListFullPath);
            writelines(fListName,logFileFullPath,"WriteMode","append");
        end

        writelines([newline,'Required program files:'],logFileFullPath,"WriteMode","append");

        for idxpList = 1:length(pList)
            pListName = pList(idxpList).Name;
            pListVer = pList(idxpList).Version;
            pListNameVer = append(pListName,"_",pListVer);
            writelines(pListNameVer,logFileFullPath,"WriteMode","append")
        end

        %% Disable diary logging
        diary off

        %% Clear figure every loop iteration; close at end of loop
        if vidNumIdx == loopEnd
            for idxFig = 1:myFigures.processNum
                idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
                close(myFigures.fig.(myFigures.handle(idxFigRef)))
            end
            myFigures.flagInLoop = false;
        else
            for idxFig = 1:myFigures.processNum
                idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
                clf(myFigures.fig.(myFigures.handle(idxFigRef)))
            end
        end

        %% Throw default error message(s)
        rethrow(ME)
    end

    % -- If no error in B2KCoordinateSystem, then: --
    % 3) Save workspace variables in .mat
    % 4) Save figure in .tif
    % 5) Store data in cell array

    %% 3) Save workspace variables for each individual run in .mat file
    wallDetFileName = append(grandchildOutputFolderNumberedName,'_wallDet');
    wallDetFileNameUnnumberedPath = fullfile(wallDetOutputFolderPath,wallDetFileName);
    wallDetFileName = B2Knextname(wallDetFileNameUnnumberedPath,'_vid_01','.mat');

    fileNamePath = fullfile(wallDetOutputFolderPath,wallDetFileName);

    % Exclude any figures or graphics objects while allowing to keep open
    % Get a list of all variables
    allvars = whos;

    % Identify the variables that ARE NOT graphics handles. This uses a regular
    % expression on the class of each variable to check if it's a graphics object
    tosave1 = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));

    % Manually remove 'myFigures' struct containing figure handles
    tosave2 = cellfun(@isempty, regexp({allvars.name},'myFigures'));

    tosavefinal = tosave1 & tosave2;

    % Pass these variable names to save
    save(fileNamePath, allvars(tosavefinal).name,'-v7.3');

    %% 4) Save figure

    for idxFig = 1:myFigures.processNum
        idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;

        % Export figure as .tif (for publication)
        fileImg = strrep(fileNamePath,'.mat','.tif');
        figNumStr = sprintf('_fig%d',idxFigRef);
        fileNameImg = insertBefore(fileImg,'.tif',figNumStr);
        % exportgraphics(myFigures.fig.(myFigures.handle(myFigures.figNum)),fileNameImg1,'Resolution',600);
        print(myFigures.fig.(myFigures.handle(idxFigRef)),fileNameImg,'-dtiffn');

        % Set CreateFcn callback
        set(myFigures.fig.(myFigures.handle(myFigures.figNum)), 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 

        % Save figure as .mat file (for reference and debugging)
        fileFig = strrep(fileNamePath,'.mat','.fig');
        fileNameFig = insertBefore(fileFig,'.fig',figNumStr);
        savefig(myFigures.fig.(myFigures.handle(idxFigRef)),fileNameFig);
    end

    %% 5) Store data in cell array
    if flagBatch == true
        cell_BatchChannelDetection{vidNumIdx,1} = myread(filteredFds.Files{vidNumIdx});
    else
        cell_BatchChannelDetection{vidNumIdx,1} = vidObj;
    end
    cell_BatchChannelDetection{vidNumIdx,2} = cell_BatchChannelDetection{vidNumIdx,1}.Name;
    cell_BatchChannelDetection{vidNumIdx,3} = tEndFnCoordSys(vidNumIdx,1);
    cell_BatchChannelDetection{vidNumIdx,4} = NaN;
    cell_BatchChannelDetection{vidNumIdx,5} = cell_BatchVideoSummarization{vidNumIdx,7};
    cell_BatchChannelDetection{vidNumIdx,6} = cell_BatchVideoSummarization{vidNumIdx,8};
    cell_BatchChannelDetection{vidNumIdx,7} = resizeScale;
    cell_BatchChannelDetection{vidNumIdx,8} = edgeMethod;
    cell_BatchChannelDetection{vidNumIdx,9} = R_highest_resolution;
    cell_BatchChannelDetection{vidNumIdx,10} = T_highest_resolution;
    cell_BatchChannelDetection{vidNumIdx,11} = H_highest_numpeaks;
    cell_BatchChannelDetection{vidNumIdx,12} = Hlines_highest_fillgap;
    cell_BatchChannelDetection{vidNumIdx,13} = Hlines_highest_minlength;
    cell_BatchChannelDetection{vidNumIdx,14} = Hlines_highest;
    cell_BatchChannelDetection{vidNumIdx,15} = Hlines_highest_num;


    %% Clear figure every loop iteration; close at end of loop
    if vidNumIdx == loopEnd
        for idxFig = 1:myFigures.processNum
            idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
            close(myFigures.fig.(myFigures.handle(idxFigRef)))
        end
        myFigures.flagInLoop = false;
    else
        for idxFig = 1:myFigures.processNum
            idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
            clf(myFigures.fig.(myFigures.handle(idxFigRef)))
        end
    end

    %% Timing end, SaveFnCoordSysData
    if flagTimeFunction == true
        tEndSaveFnCoordSysData(vidNumIdx,1) = toc(tStartSaveFnCoordSysData);
    end

    cell_BatchChannelDetection{vidNumIdx,4} = tEndSaveFnCoordSysData(vidNumIdx,1);

    %% B2KParticleTracking - Detect and track particle

    % Plot variables
    tile_x_pos = -200;
    tile_y_pos = 30;

    %% Timing start, FnParticleTrack
    if flagTimeFunction == true
        tStartFnParticleTrack = tic;
    end

    try
        [myFigures,time_ms,position_x_CS,velocity_x_CS,acceleration_x_CS,position_y_CS,velocity_y_CS,acceleration_y_CS,meanThroughput_x_CS] = B2KParticleTrackingBGM(vidNumIdx,vidObj,myFigures,onePTrackOutputFolderPath,L_start,R_end,tFrames,nFrames,I,CS,x_pos_caption,y_pos_caption,EXP_actualframerate,flagDebug,flagRunWithDisplay);
        
        %% Timing end, FnParticleTrack
        if flagTimeFunction == true
            tEndFnParticleTrack(vidNumIdx,1) = toc(tStartFnParticleTrack);
        end

        %% Timing start, SaveFnParticleTrackData
        if flagTimeFunction == true
            tStartSaveFnParticleTrackData = tic;
        end
    catch ME
        %% Timing end, FnParticleTrack
        if flagTimeFunction == true
            tEndFnParticleTrack(vidNumIdx,1) = toc(tStartFnParticleTrack);
        end

        %% Timing start, SaveFnParticleTrackData
        if flagTimeFunction == true
            tStartSaveFnParticleTrackData = tic;
        end

        %% Manually add appropriate error message(s) from function call(s)
        % Construct a detailed error message including the stack trace
        errorMsg = sprintf('\n%s\n\n', ME.message);
    
        for k = 1:length(ME.stack)-1
            % Read the file and get the specific line where the error occurred
            fileContent = fileread(ME.stack(k).file);
            fileLines = strsplit(fileContent, '\n');
            errorLine = fileLines{ME.stack(k).line};
    
            % Create clickable bold link using MATLAB's hyperlink format with HTML bold tags
            errorMsg = sprintf('%sError in <a href="matlab: opentoline(''%s'',%d)">%s</a> (line %d):\n%s\n', ...
                               errorMsg, ME.stack(k).file, ME.stack(k).line, ME.stack(k).name, ME.stack(k).line, errorLine);
        end
    
        % Use the error function to display the error message in red and stop execution
        error(errorMsg);

        %%
        % -- If error in B2KParticleTracking, then: --
        % 1) Save table results to .mat
        % 2) Save used functions to log .txt

        table_BatchParticleTracking = cell2table(cell_BatchParticleTracking,...
            "VariableNames",...
            ["Video Reader Object" ...                      % 1
            "Video File Name" ...                           % 2
            "Time Particle Tracking (s)"                    % 3
            "Time Save Particle Tracking (s)"               % 4
            "Start Frame" ...                               % 5
            "End Frame" ...                                 % 6
            "Time Array (s)" ...                            % 7
            "x Position Array (mm)" ...                     % 8
            "x Velocity Array (m/s)" ...                    % 9
            "x Acceleration Array (m/s^2)" ...              % 10
            "y Position Array (m)" ...                      % 11
            "y Velocity Array (m/s)" ...                    % 12
            "y Acceleration Array (m/s^2)" ...              % 13
            "Mean Particle Throughput (m/s)"]);             % 14

        %% 1) Save table results to .mat
        trackTableOutputName = append("table_",grandchildOutputFolderNumberedName,"_Track",".mat");
        trackTableOutputPath = fullfile(wallDetOutputFolderPath,trackTableOutputName);

        % save(fileNamePath);
        save(trackTableOutputPath,"table_BatchParticleTracking",'-v7.3');

        %% 2) Save used functions to log .txt file

        currScriptNameWithExt = append(currScriptName,'.m');
        [fList,pList] = matlab.codetools.requiredFilesAndProducts(currScriptNameWithExt);

        logFileName = append(grandchildOutputFolderNumberedName,'_log','.txt');
        logFileFullPath = fullfile(grandchildOutputFolderPath,logFileName);
        writelines([newline,'Function run (function_YYMMDD_instance):'],logFileFullPath);
        writelines(grandchildOutputFolderNumberedName,logFileFullPath,"WriteMode","append");
        
        writelines([newline,'Main function:'],logFileFullPath,"WriteMode","append");
        writelines(currScriptNameWithExt,logFileFullPath,"WriteMode","append");

        writelines([newline,'Settings:'],logFileFullPath,"WriteMode","append");
        writelines(sprintf('Added frames: %d',userOption.addedFrames),logFileFullPath,"WriteMode","append");
        writelines(sprintf('Video number start: %d',loopStart),logFileFullPath,"WriteMode","append");
        writelines(sprintf('Video number end: %d',loopEnd),logFileFullPath,"WriteMode","append");

        writelines([newline,'Flags: '],logFileFullPath,"WriteMode","append");
        writelines(['Drive priority: ',userOption.drivePriority],logFileFullPath,"WriteMode","append");
        writelines(['Processing mode: ',userOption.processingMode],logFileFullPath,"WriteMode","append");
        writelines(['Run mode: ',userOption.runMode],logFileFullPath,"WriteMode","append");
        writelines(['Display mode: ',userOption.displayMode],logFileFullPath,"WriteMode","append");
        writelines(['Save location: ',userOption.saveLocation],logFileFullPath,"WriteMode","append");
        writelines(['Time process function: ',userOption.evalTimeFunction],logFileFullPath,"WriteMode","append");
        writelines(['Time saving process function data: ',userOption.evalTimeSaveFunctionData],logFileFullPath,"WriteMode","append");
        writelines(['Data path override: ',userOption.dataPathOverride],logFileFullPath,"WriteMode","append");
        writelines(['Drive override: ',userOption.driveOverride],logFileFullPath,"WriteMode","append");
        writelines(['Single video override: ',userOption.singleVidOverride],logFileFullPath,"WriteMode","append");

        writelines([newline,'Functions used:'],logFileFullPath,"WriteMode","append");
        for idxfList = 1:length(fList)
            fListFullPath = fList{idxfList};
            [fListFilePath,fListName,FListExt] = fileparts(fListFullPath);
            writelines(fListName,logFileFullPath,"WriteMode","append");
        end

        writelines([newline,'Required program files:'],logFileFullPath,"WriteMode","append");

        for idxpList = 1:length(pList)
            pListName = pList(idxpList).Name;
            pListVer = pList(idxpList).Version;
            pListNameVer = append(pListName,"_",pListVer);
            writelines(pListNameVer,logFileFullPath,"WriteMode","append")
        end

        %% Disable diary logging
        diary off

        %% Clear figure every loop iteration; close at end of loop
        if vidNumIdx == loopEnd
            for idxFig = 1:myFigures.processNum
                idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
                close(myFigures.fig.(myFigures.handle(idxFigRef)))
            end
            myFigures.flagInLoop = false;
        else
            for idxFig = 1:myFigures.processNum
                idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
                clf(myFigures.fig.(myFigures.handle(idxFigRef)))
            end
        end

        %% Throw default error message(s)
        rethrow(ME);
    end

    %% 3) Save workspace variables for each individual run in .mat file
    onePTrackFileName = append(grandchildOutputFolderNumberedName,'_1pTracking');
    onePTrackFileNameUnnumberedPath = fullfile(onePTrackOutputFolderPath,onePTrackFileName);
    onePTrackFileName = B2Knextname(onePTrackFileNameUnnumberedPath,'_vid_01','.mat');

    fileNamePath = fullfile(onePTrackOutputFolderPath,onePTrackFileName);

    % Exclude any figures or graphics objects while allowing to keep open
    % Get a list of all variables
    allvars = whos;

    % Identify the variables that ARE NOT graphics handles. This uses a regular
    % expression on the class of each variable to check if it's a graphics object
    tosave1 = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));

    % Manually remove 'myFigures' struct containing figure handles
    tosave2 = cellfun(@isempty, regexp({allvars.name},'myFigures'));

    tosavefinal = tosave1 & tosave2;

    % Pass these variable names to save
    save(fileNamePath, allvars(tosavefinal).name,'-v7.3');

    %% 4) Save figure

    for idxFig = 1:myFigures.processNum
        idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;

        % Export figure as .tif (for publication)
        fileImg = strrep(fileNamePath,'.mat','.tif');
        figNumStr = sprintf('_fig%d',idxFigRef);
        fileNameImg = insertBefore(fileImg,'.tif',figNumStr);
        % exportgraphics(myFigures.fig.(myFigures.handle(myFigures.figNum)),fileNameImg1,'Resolution',600);
        print(myFigures.fig.(myFigures.handle(idxFigRef)),fileNameImg,'-dtiffn');

        % Set CreateFcn callback
        set(myFigures.fig.(myFigures.handle(myFigures.figNum)), 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 

        % Save figure as .mat file (for reference and debugging)
        fileFig = strrep(fileNamePath,'.mat','.fig');
        fileNameFig = insertBefore(fileFig,'.fig',figNumStr);
        savefig(myFigures.fig.(myFigures.handle(idxFigRef)),fileNameFig);
    end

    %% 5) Store data in cell array
    cell_BatchParticleTracking{vidNumIdx,1} = vidObj;
    cell_BatchParticleTracking{vidNumIdx,2} = vidObj.Name;
    cell_BatchParticleTracking{vidNumIdx,3} = tEndFnParticleTrack(vidNumIdx,1);
    cell_BatchParticleTracking{vidNumIdx,4} = NaN; %Data stored later
    cell_BatchParticleTracking{vidNumIdx,5} = cell_BatchVideoSummarization{vidNumIdx,7};
    cell_BatchParticleTracking{vidNumIdx,6} = cell_BatchVideoSummarization{vidNumIdx,8};
    cell_BatchParticleTracking{vidNumIdx,7} = time_ms;
    cell_BatchParticleTracking{vidNumIdx,8} = position_x_CS;
    cell_BatchParticleTracking{vidNumIdx,9} = velocity_x_CS;
    cell_BatchParticleTracking{vidNumIdx,10} = acceleration_x_CS;
    cell_BatchParticleTracking{vidNumIdx,11} = position_y_CS;
    cell_BatchParticleTracking{vidNumIdx,12} = velocity_y_CS;
    cell_BatchParticleTracking{vidNumIdx,13} = acceleration_y_CS;
    cell_BatchParticleTracking{vidNumIdx,14} = meanThroughput_x_CS;

    %% Clear figure every loop iteration; close at end of loop
    if vidNumIdx == loopEnd
        for idxFig = 1:myFigures.processNum
            idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
            close(myFigures.fig.(myFigures.handle(idxFigRef)))
        end
        myFigures.flagInLoop = false;
    else
        for idxFig = 1:myFigures.processNum
            idxFigRef = myFigures.figNumArray(end) - myFigures.processNum + idxFig;
            clf(myFigures.fig.(myFigures.handle(idxFigRef)))
        end

        %% Reset figure counter to initial figure number for loop
        myFigures.figNum = myFigures.initFigNumLoop;
    end

    %% Timing end, SaveFnParticleTrackData
    if flagTimeFunction == true
        tEndSaveFnParticleTrackData(vidNumIdx,1) = toc(tStartSaveFnParticleTrackData);
    end

    cell_BatchParticleTracking{vidNumIdx,4} = tEndSaveFnParticleTrackData(vidNumIdx,1);
end

%% Video summarization - log and data table
table_BatchVideoSummarization = cell2table(cell_BatchVideoSummarization,...
    "VariableNames",...
    ["Video Reader Object" ...
    "Video File Name" ...
    "Time Video Summary (s)" ...
    "Time Save Video Summary (s)" ...
    "Error Count" ...
    "Pass Count" ...
    "Start Frame" ...
    "End Frame" ...
    "Number of Frames to Analyze" ...
    "Total Number of Frames" ...
    "Flag p-Chip Present" ...
    "Flag p-Chip Stuck" ...
    "Flag no p-Chip" ...
    "Flag Video Corrupt"]);

%% Video summarization - save table results to .mat file
vidSumTableOutputName = append("table_",vidSumFileName,".mat");
vidSumTableOutputPath = fullfile(vidSumOutputFolderPath,vidSumTableOutputName);

% save(fileNamePath);
save(vidSumTableOutputPath,"table_BatchVideoSummarization",'-v7.3');

%% Channel detection - log and data table
table_BatchChannelDetection = cell2table(cell_BatchChannelDetection,...
    "VariableNames",...
    ["Video Reader Object" ...
    "Video File Name" ...
    "Time Coordinate System (s)" ...
    "Time Save Coordinate System (s)" ...
    "Start Frame" ...
    "End Frame" ...
    "Resize Scale" ...
    "Edge Method" ...
    "Hough Rho Resolution" ...
    "Hough Theta Resolution" ...
    "Hough Number of Peaks" ...
    "Houghlines Fill Gap" ...
    "Houghlines Minimum Length" ...
    "Houghlines" ...
    "Houghlines Number Detected"]);

%% Channel detection - save table results to .mat file
wallDetTableOutputName = append("table_",wallDetFileName,".mat");
wallDetTableOutputPath = fullfile(wallDetOutputFolderPath,wallDetTableOutputName);

% save(fileNamePath);
save(wallDetTableOutputPath,"table_BatchChannelDetection",'-v7.3');

%% Particle tracking - log and data table
table_BatchParticleTracking = cell2table(cell_BatchParticleTracking,...
    "VariableNames",...
    ["Video Reader Object" ...                      % 1
    "Video File Name" ...                           % 2
    "Time Particle Tracking (s)" ...                % 3
    "Time Save Particle Tracking (s)" ...           % 4
    "Start Frame" ...                               % 5
    "End Frame" ...                                 % 6
    "Time Array (s)" ...                            % 7
    "x Position Array (mm)" ...                     % 8
    "x Velocity Array (m/s)" ...                    % 9
    "x Acceleration Array (m/s^2)" ...              % 10
    "y Position Array (mm)" ...                     % 11
    "y Velocity Array (m/s)" ...                    % 12
    "y Acceleration Array (m/s^2)" ...              % 13
    "Mean Particle Throughput (m/s)"]);             % 14

%% Particle tracking
trackTableOutputName = append("table_",grandchildOutputFolderNumberedName,"_Track",".mat");
trackTableOutputPath = fullfile(onePTrackOutputFolderPath,trackTableOutputName);

% save(fileNamePath);
save(trackTableOutputPath,"table_BatchParticleTracking",'-v7.3');

%% Save used functions to log .txt file
currScriptNameWithExt = append(currScriptName,'.m');
[fList,pList] = matlab.codetools.requiredFilesAndProducts(currScriptNameWithExt);

logFileName = append(grandchildOutputFolderNumberedName,'_log','.txt');
logFileFullPath = fullfile(grandchildOutputFolderPath,logFileName);
writelines([newline,'Function run (function_YYMMDD_instance):'],logFileFullPath);
writelines(grandchildOutputFolderNumberedName,logFileFullPath,"WriteMode","append");

writelines([newline,'Main function:'],logFileFullPath,"WriteMode","append");
writelines(currScriptNameWithExt,logFileFullPath,"WriteMode","append");

writelines([newline,'Settings:'],logFileFullPath,"WriteMode","append");
writelines(sprintf('Added frames: %d',userOption.addedFrames),logFileFullPath,"WriteMode","append");
writelines(sprintf('Video number start: %d',loopStart),logFileFullPath,"WriteMode","append");
writelines(sprintf('Video number end: %d',loopEnd),logFileFullPath,"WriteMode","append");

writelines([newline,'Flags: '],logFileFullPath,"WriteMode","append");
writelines(['Drive priority: ',userOption.drivePriority],logFileFullPath,"WriteMode","append");
writelines(['Processing mode: ',userOption.processingMode],logFileFullPath,"WriteMode","append");
writelines(['Run mode: ',userOption.runMode],logFileFullPath,"WriteMode","append");
writelines(['Display mode: ',userOption.displayMode],logFileFullPath,"WriteMode","append");
writelines(['Save location: ',userOption.saveLocation],logFileFullPath,"WriteMode","append");
writelines(['Time process function: ',userOption.evalTimeFunction],logFileFullPath,"WriteMode","append");
writelines(['Time saving process function data: ',userOption.evalTimeSaveFunctionData],logFileFullPath,"WriteMode","append");
writelines(['Data path override: ',userOption.dataPathOverride],logFileFullPath,"WriteMode","append");
writelines(['Drive override: ',userOption.driveOverride],logFileFullPath,"WriteMode","append");
writelines(['Single video override: ',userOption.singleVidOverride],logFileFullPath,"WriteMode","append");

writelines([newline,'Functions used:'],logFileFullPath,"WriteMode","append");
for idxfList = 1:length(fList)
    fListFullPath = fList{idxfList};
    [fListFilePath,fListName,FListExt] = fileparts(fListFullPath);
    writelines(fListName,logFileFullPath,"WriteMode","append");
end

writelines([newline,'Required program files:'],logFileFullPath,"WriteMode","append");

for idxpList = 1:length(pList)
    pListName = pList(idxpList).Name;
    pListVer = pList(idxpList).Version;
    pListNameVer = append(pListName,"_",pListVer);
    writelines(pListNameVer,logFileFullPath,"WriteMode","append")
end

%% Clear functions
clear functions

%% Timing - whole script, end
tScriptEnd = toc(tScriptStart) % pair 1: toc

%% Disable diary logging
diary off

%% Profile code [end]

% profile viewer
% profsave

% Stops profiler and displays struct
profileHandle = profile('info'); 

% Specify full absolute path to save profiler .MAT
profileOutputName = append(grandchildOutputFolderNumberedName,'_profile','.mat');
profileOutputFullPath = fullfile(grandchildOutputFolderPath,profileOutputName);

% Save profile data
save(profileOutputFullPath,"profileHandle",'-v7.3')
clear profileHandle

% To check Profiler results, 'load' and use 'profview' to load profile handle, Example:
% load("F:\OutputFolder\Batch1pTracking\Batch1pTracking_230816_07\Batch1pTracking_230816_07_profile.mat")
% profview(0,profileHandle)

%% ===== LOCAL FUNCTIONS =====

%% Function: vtest - custom reader
function vtest = myread(file)
vtest = VideoReader(file);
end

%% -------- Additional local helper functions --------
function root = getProjectRootOrHere()
% Try MATLAB Project root; fall back to this file’s folder
try
    p = matlab.project.currentProject();
    root = p.RootFolder;
    if ~isfolder(root), error('Invalid project root.'); end
catch
    here = fileparts(currentFileFullpath());
    if isempty(here), here = pwd; end
    root = here;
end
end

function out = getDownloadsFolder()
% Cross-platform "Downloads" folder path with simple fallbacks
if ispc
    home = getenv('USERPROFILE');
    if isempty(home), home = char(java.lang.System.getProperty('user.home')); end
    out = fullfile(home,'Downloads');
elseif ismac || isunix
    home = getenv('HOME');
    if isempty(home), home = char(java.lang.System.getProperty('user.home')); end
    out = fullfile(home,'Downloads');
else
    out = fullfile(pwd,'Downloads'); % last-resort fallback
end
mkIfMissing(out);
end

function mkIfMissing(p)
if exist(p,'dir') ~= 7
    mkdir(p);
end
end

function name = currentScriptNameOr(defaultName)
% Robust name for scripts/functions (mfilename can be empty for scripts/livescripts)
name = mfilename;
if isempty(name)
    [~, name, ~] = fileparts(currentFileFullpath());
    if isempty(name), name = defaultName; end
end
end

function f = currentFileFullpath()
% Full path of current file (works for scripts/functions)
f = '';
st = dbstack('-completenames');
if numel(st) >= 1
    f = st(end).file;  % bottom of stack is the script being run
elseif usejava('desktop')
    try
        doc = matlab.desktop.editor.getActive;
        if ~isempty(doc) && doc.Opened
            f = doc.Filename;
        end
    catch
        % leave empty
    end
end
end

function paths = ensureCharCell(paths)
% Ensure input is a cell array of char vectors (for older fileDatastore versions)
if isstring(paths)
    paths = cellstr(paths);              % string array -> cellstr
elseif ischar(paths)
    paths = {paths};                     % char -> {char}
elseif iscell(paths)
    % If it's a cell of strings, convert to cell of char
    if ~isempty(paths) && isstring(paths{1})
        paths = cellfun(@char, paths, 'UniformOutput', false);
    end
else
    error('Unsupported path type for fileDatastore.');
end
end
%% Batch1pTracking04_orig_Final.m - [Script] Batch single p-Chip tracking
%
% Author: Raymond Yeung
% Release date: 2023
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
% profile on -historysize 10000000
% profile on -historysize 20000000 [Output: Flame graph is not available because the number of function calls exceeds the current profiler history size of 20000000.]
profile on -historysize 100000000

%% Timing - whole script, start
tScriptStart = tic; % pair 1: tic

%% Choose User
userOption.User = 'UserOne';

%% User Options (storing all variables that user may change to tweak code)
if strcmp(userOption.User,'UserOne')
    % [SETTINGS: Program modes]
    % Drive priority: external, internal (legacy mode only)
    % Processing modes = batch, single
    % Run modes: debug, final
    % Display modes: runWithDisplay, runWithoutDisplay
    % Save location: defaultVideoLocation, anotherLocation
    userOption.drivePriority = 'internal';
    userOption.processingMode = 'batch';
    userOption.runMode = 'debug';
    userOption.displayMode = 'runWithDisplay';
    userOption.saveLocation = 'defaultVideoLocation';

    % [NEW] Preferred data source: 'project' (recommended) or 'drives' (legacy)
    userOption.dataSourceMode = 'project';

    % [SETTINGS (OPTIONAL): Stored Data Path] — keep empty for packaging
    userOption.dataPathOverride = '';

    % [SETTINGS (OPTIONAL): Drive]  (legacy search; keep empty for portability)
    userOption.driveOverride = '';

    % [SETTINGS (OPTIONAL): Single video analysis]
    userOption.singleVidOverride = '';

    % [SETTINGS: Video Summarization]
    userOption.addedFrames = 25;

else
    error('Unknown user name.');
end

%% Figure number definition using structures containing fields with dynamic expression
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
[flagInternalDrive,flagBatch,flagDebug,flagRunWithDisplay,flagSaveLocDefault,flagDataPathOverride,flagDriveOverride,flagSingleVidOverride] = ...
    B2KFlags(userOption.drivePriority,userOption.processingMode,userOption.runMode,userOption.displayMode,userOption.saveLocation,userOption.dataPathOverride,userOption.driveOverride,userOption.singleVidOverride);

% Local flag for project-root data source (for packaged project)
flagProjectSource = strcmpi(userOption.dataSourceMode,'project');

%% Find folder storing data videos and set up output folder

% Saved default drives / data videos folder / output folder locations
ExpInternalDriveNames = ["INT_DRIVE_A","INT_DRIVE_B","INT_DRIVE_C"];
ExpExternalDriveNames = ["EXT_DRIVE_A","EXT_DRIVE_B","EXT_DRIVE_C"];

% List all physical drives on computer (used only in legacy mode)
infoTable = B2KlistPhysicalDrives;

% Spit warning if user provided options for both the dataPathOverride and driveOverride
if flagDataPathOverride && flagDriveOverride
    warning('Both dataPathOverride and driveOverride provided. Defaulting to dataPathOverride ...')
end

if flagDataPathOverride
    % 1) Hard override of data path takes precedence
    selpath = userOption.dataPathOverride;
    fprintf('Analyzing data from chosen (overridden) location, %s ...\n', selpath);

elseif flagProjectSource
    % 2) Project-root anchored (for packaged projects)
    projRoot = getProjectRootOrHere();                 % absolute path to the MATLAB Project root (or this file’s folder)
    dv = fullfile(projRoot,'Data_Videos');             % videos expected in <project root>/Data_Videos/*.avi
    if exist(dv,'dir') ~= 7
        warning('[B2K] Expected "<project>/Data_Videos" not found. Please select the working folder that contains "Data_Videos".')
        selpath = uigetdir;
        if isequal(selpath,0), error('No folder selected.'); end
        fprintf('Analyzing data from chosen (GUI) location, %s ...\n', selpath);
    else
        selpath = projRoot;                             % anchor to project root for reading videos
        fprintf('Analyzing data from Project root: %s\n', selpath);
    end

else
    % 3) Legacy drive search
    if ~flagDriveOverride
        if flagInternalDrive
            ExpDriveNames = [ExpInternalDriveNames,ExpExternalDriveNames];
        else
            ExpDriveNames = [ExpExternalDriveNames,ExpInternalDriveNames];
        end

        matchDriveName = [];
        for d = 1:numel(ExpDriveNames)
            if ispc
                matchDriveName = infoTable.DeviceID(strcmp(infoTable.VolumeName,ExpDriveNames(d)),:);
            elseif ismac
                matchDriveName = fullfile('/Volumes',infoTable.VolumeName(strcmp(infoTable.VolumeName,ExpDriveNames(d)),:));
            else
                error('Have not tested in Linux')
            end
            if ~isempty(matchDriveName), break; end
        end

        if ~isempty(matchDriveName)
            selpath = matchDriveName;
            if flagInternalDrive
                if d <= numel(ExpInternalDriveNames)
                    fprintf('Analyzing data from internal drive, %s%s ...\n', matchDriveName,ExpDriveNames(d));
                else
                    fprintf('Analyzing data from external drive, %s%s ...\n', matchDriveName,ExpDriveNames(d));
                end
            else
                if d > numel(ExpInternalDriveNames)
                    fprintf('Analyzing data from internal drive, %s%s ...\n', matchDriveName,ExpDriveNames(d));
                else
                    fprintf('Analyzing data from external drive, %s%s ...\n', matchDriveName,ExpDriveNames(d));
                end
            end
        else
            warning('[B2K] Default drives not found. Please select the working folder. The working folder should have a subfolder titled "Data_Videos" containing the experimental videos.')
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

%5a Current script name (robust)
currScriptName = currentScriptNameOr('MainScript');

%5b Date ('yyMMdd')
currDate = char(datetime('now','Format','yyMMdd'));

%5 Combine current script name and date
currScriptNameAddDate = append(currScriptName,'_',currDate);

%4 Child output folder (extracted from script name)
childOutputFolderName = char(regexp(currScriptName,'[a-zA-Z]+\w*[a-zA-Z]+','match','once'));
if isempty(childOutputFolderName), childOutputFolderName = currScriptName; end

%3 Parent output folder
parentOutputFolderName = 'OutputFolder';

%4 Paths
childOutputFolderPath = fullfile(outputBase,parentOutputFolderName,childOutputFolderName);
grandchildOutputFolderUnnumberedPath = fullfile(childOutputFolderPath,currScriptNameAddDate);
grandchildOutputFolderNumberedName   = B2Knextname(grandchildOutputFolderUnnumberedPath,'_01','');
grandchildOutputFolderPath           = fullfile(outputBase,parentOutputFolderName,childOutputFolderName,grandchildOutputFolderNumberedName);

% Ensure folders exist
mkIfMissing(childOutputFolderPath);
mkIfMissing(grandchildOutputFolderPath);

% Sub-output folders
vidSumOutputFolderPath  = fullfile(grandchildOutputFolderPath,'VideoSummarization_Output');   mkIfMissing(vidSumOutputFolderPath);
wallDetOutputFolderPath = fullfile(grandchildOutputFolderPath,'ChannelWallDetection_Output'); mkIfMissing(wallDetOutputFolderPath);
onePTrackOutputFolderPath = fullfile(grandchildOutputFolderPath,'OnePTracking_Output');       mkIfMissing(onePTrackOutputFolderPath);

%% Create diary file for logging commands, keyboard inputs, text output
diaryOutputName       = append(grandchildOutputFolderNumberedName,'_diary');
diaryOutputFolderPath = fullfile(grandchildOutputFolderPath,diaryOutputName);
diary(diaryOutputFolderPath)

%% GUI functions
% uigetfile, uiputfile, uigetdir, uiopen, uisave

%% Video datastore: Single or Batch
% Project mode -> use ONLY "<project root>/Data_Videos"
if flagProjectSource
    baseVideos = fullfile(selpath,'Data_Videos');   % <project root>/Data_Videos
    if exist(baseVideos,'dir') ~= 7
        error('Project mode: folder "%s" not found. Expected at "<project root>/Data_Videos".', baseVideos);
    end

    if flagBatch
        fds = fileDatastore(baseVideos,'ReadFcn',@myread,'IncludeSubfolders',true,'FileExtensions','.avi');
        vidPathStack = fds.Files;
        if isempty(vidPathStack)
            error('No .avi videos found under "%s".', baseVideos);
        end
    else
        if flagSingleVidOverride && ~isempty(userOption.singleVidOverride)
            warning('Single video mode: using provided singleVidOverride ...')
            vidPathStack = userOption.singleVidOverride;
        else
            warning('Single video mode: please select a video from the "Data_Videos" folder.')
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

            uialert(figSingleVidAlert,'Program currently set to analyze a single video. Please select a video from the "Data_Videos" folder.','Single video processing mode','Modal',false,'Icon','info','CloseFcn','uiresume(figSingleVidAlert)');
            uiwait(figSingleVidAlert);
            close(figSingleVidAlert);

            [vidName,vidPath] = uigetfile({'*.avi','AVI Video (*.avi)'},'Select video',baseVideos);
            if isequal(vidName,0), error('No video selected.'); end
            vidPathStack = fullfile(vidPath,vidName);
        end
    end

else
    % Legacy mode: look under <selpath>/Data_Videos (single folder)
    loc = fullfile(selpath,'Data_Videos');
    if flagBatch
        fds = fileDatastore(loc,'ReadFcn',@myread,'IncludeSubfolders',true,'FileExtensions','.avi');
        vidPathStack = fds.Files;
        if isempty(vidPathStack)
            error('No .avi videos found under "%s".', loc);
        end
    else
        if flagSingleVidOverride && ~isempty(userOption.singleVidOverride)
            warning('Program currently set to analyze a single video. Using singleVidOverride ...')
            vidPathStack = userOption.singleVidOverride;
        else
            warning('Program currently set to analyze a single video. Please select a video from the "Data_Videos" folder.')
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

            uialert(figSingleVidAlert,'Program currently set to analyze a single video. Please select a video from the "Data_Videos" folder.','Single video processing mode','Modal',false,'Icon','info','CloseFcn','uiresume(figSingleVidAlert)');
            uiwait(figSingleVidAlert);
            close(figSingleVidAlert);

            [vidName,vidPath] = uigetfile({'*.avi','AVI Video (*.avi)'},'Select video',loc);
            if isequal(vidName,0), error('No video selected.'); end
            vidPathStack = fullfile(vidPath,vidName);
        end
    end
end

%% Single/batch data table setup

DataTable = table(); %STILL NEED TO ADD ALL DATA TO THIS DATATABLE

%% Video Summarization setup
addedFrames = userOption.addedFrames;

cell_BatchVideoSummarization = cell(size(vidPathStack,1),31);
vidSumPassCount = 0;
vidSumErrCount = 0;

%% Channel Wall Detection setup

idx_cell_BatchChannelDetection = 13;
cell_BatchChannelDetection = cell(size(vidPathStack,1),idx_cell_BatchChannelDetection);

cell_AllData = cell(2,2);

HLcorrectionidx = 0;

%% Particle Tracking setup

idx_cell_ParticleTracking = 13;
cell_ParticleTracking = cell(size(vidPathStack,1),idx_cell_ParticleTracking);

%% Loop

loopStart = 1;
loopEnd = 1;
% loopEnd = length(fds.Files)

%% [FIGURE TEST]
myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);

% Set CreateFcn callback
set(myFigures.fig.(myFigures.handle(myFigures.figNum)), 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 

%
clear B2KMyFigures
% save figure
close(myFigures.fig.(myFigures.handle(myFigures.figNum)))

for i = loopStart:loopEnd %change length(fds.Files) to accomodate single and batch
    %% Update flag loop status
    myFigures.flagInLoop = true;

    if i == loopStart
        myFigures.initFigNumLoop = myFigures.figNum;
    end

    %% Read video to create video object and vid name
    if flagBatch == true
        vidPath = fds.Files{i};
        vidObj = myread(vidPath);
        [vidFolder,vidNameNoExt,ext] = fileparts(vidPath);
        vid = append(vidNameNoExt,ext);
    else %if flagBatch == false
        vidPath = vidPathStack;
        vidObj = VideoReader(vidPath);
        [vidFolder,vidNameNoExt,ext] = fileparts(vidPath);
        vid = append(vidNameNoExt,ext);
    end

    %% Number
    nFrames = vidObj.NumFrames;

    %% B2KVideoSummary - Video summarization to identify start and end frames to analyze

    try
        [L_start,R_end,tFrames,nFrames,...
            flag_pChip_present,flag_pChip_stuck,flag_no_pChip,flag_video_corrupt,...
            myFigures,...
            thresh_median_idx, thresh_median_pks, thresh_median_count, thresh_mean_idx, thresh_mean_pks, thresh_mean_count,...
            thresh_quartiles_idx, thresh_quartiles_pks, thresh_quartiles_count, thresh_gesd_idx, thresh_gesd_pks, thresh_gesd_count,...
            TF,S1,S2,...
            ipt,residual,...
            L_startBlockStd,R_endBlockStd] = B2KVideoSummary01(vidObj,myFigures,flagDebug,flagRunWithDisplay,addedFrames);
    catch ME
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
            ["Video Reader Object" ...
            "Video File Name" ...
            "Error Count" ...
            "Pass Count" ...
            "Start Frame" ...
            "End Frame" ...
            "Number of Frames to Analyze" ...
            "Total Number of Frames" ...
            "Flag p-Chip Present" ...
            "Flag p-Chip Stuck" ...
            "Flag no p-Chip" ...
            "Flag Video Corrupt" ...
            "Threshold Median Index" ...
            "Threshold Median Peaks" ...
            "Threshold Median Count" ...
            "Threshold Mean Index" ...
            "Threshold Mean Peaks" ...
            "Threshold Mean Count" ...
            "Threshold Quartiles Index" ...
            "Threshold Quartiles Peaks" ...
            "Threshold Quartiles Count" ...
            "Threshold GESD Index" ...
            "Threshold GESD Peaks" ...
            "Threshold GESD Count" ...
            "TF" ...
            "S1" ...
            "S2" ...
            "ipt" ...
            "residual" ...
            "Start StDev" ...
            "End StDev"]);

        %% 1) Save table results to .mat
        vidSumTableOutputName = append("table_",grandchildOutputFolderNumberedName,"_vidSum",".mat");
        vidSumTableOutputPath = fullfile(vidSumOutputFolderPath,vidSumTableOutputName);

        save(vidSumTableOutputPath,"table_BatchVideoSummarization");

        %% 2) Save used functions to log .txt file

        currScriptNameWithExt = append(currScriptName,'.m');
        [fList,pList] = matlab.codetools.requiredFilesAndProducts(currScriptNameWithExt);

        logFileName = append(grandchildOutputFolderNumberedName,'_log','.txt');
        logFileFullPath = fullfile(vidSumOutputFolderPath,logFileName);
        writelines(grandchildOutputFolderNumberedName,logFileFullPath);

        writelines(currScriptNameWithExt,logFileFullPath,"WriteMode","append");
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
        if i == loopEnd
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
    save(fileNamePath, allvars(tosavefinal).name);

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
    cell_BatchVideoSummarization{i,1} = vidObj;
    cell_BatchVideoSummarization{i,2} = vidObj.Name;
    cell_BatchVideoSummarization{i,3} = NaN;
    cell_BatchVideoSummarization{i,4} = NaN;
    cell_BatchVideoSummarization{i,5} = L_start;
    cell_BatchVideoSummarization{i,6} = R_end;
    cell_BatchVideoSummarization{i,7} = tFrames;
    cell_BatchVideoSummarization{i,8} = nFrames;
    cell_BatchVideoSummarization{i,9} = flag_pChip_present;
    cell_BatchVideoSummarization{i,10} = flag_pChip_stuck;
    cell_BatchVideoSummarization{i,11} = flag_no_pChip;
    cell_BatchVideoSummarization{i,12} = flag_video_corrupt;
    cell_BatchVideoSummarization{i,13} = thresh_median_idx;
    cell_BatchVideoSummarization{i,14} = thresh_median_pks;
    cell_BatchVideoSummarization{i,15} = thresh_median_count;
    cell_BatchVideoSummarization{i,16} = thresh_mean_idx;
    cell_BatchVideoSummarization{i,17} = thresh_mean_pks;
    cell_BatchVideoSummarization{i,18} = thresh_mean_count;
    cell_BatchVideoSummarization{i,19} = thresh_quartiles_idx;
    cell_BatchVideoSummarization{i,20} = thresh_quartiles_pks;
    cell_BatchVideoSummarization{i,21} = thresh_quartiles_count;
    cell_BatchVideoSummarization{i,22} = thresh_gesd_idx;
    cell_BatchVideoSummarization{i,23} = thresh_gesd_pks;
    cell_BatchVideoSummarization{i,24} = thresh_gesd_count;
    cell_BatchVideoSummarization{i,25} = TF;
    cell_BatchVideoSummarization{i,26} = S1;
    cell_BatchVideoSummarization{i,27} = S2;
    cell_BatchVideoSummarization{i,28} = ipt;
    cell_BatchVideoSummarization{i,29} = residual;
    cell_BatchVideoSummarization{i,30} = L_startBlockStd;
    cell_BatchVideoSummarization{i,31} = R_endBlockStd;

    %%
    if isnan(flag_pChip_present) || ~isnan(flag_pChip_stuck) || ~isnan(flag_no_pChip) || ~isnan(flag_video_corrupt)
        vidSumErrCount = vidSumErrCount + 1;
        cell_BatchVideoSummarization{i,3} = vidSumErrCount;
        %i %report index to keep track of computing time
        i
        clf(myFigures.fig.(myFigures.handle(myFigures.figNum)))
        %         close(fig1)

        % if i ~= loopEnd
        %     continue
        % else

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
        save(fileNamePath, allvars(tosavefinal).name);

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
        cell_BatchChannelDetection{i,1} = myread(fds.Files{i});
        cell_BatchChannelDetection{i,2} = cell_BatchChannelDetection{i,1}.Name;

        for j = 3:idx_cell_BatchChannelDetection
            cell_BatchChannelDetection{i,j} = NaN;
        end

        %% Clear figure every loop iteration; close at end of loop
        if i == loopEnd
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

        continue

    end

    %%
    vidSumPassCount = vidSumPassCount + 1;
    cell_BatchVideoSummarization{i,4} = vidSumPassCount;

    %% Clear figure every loop iteration; close at end of loop
    if i == loopEnd
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

    %% Read results from BatchVideoSummarization and skip videos with errors
    if isnan(cell_BatchVideoSummarization{i,4})

        cell_BatchChannelDetection{i,1} = myread(fds.Files{i});
        cell_BatchChannelDetection{i,2} = cell_BatchChannelDetection{i,1}.Name;

        for j = 3:idx_cell_BatchChannelDetection
            cell_BatchChannelDetection{i,j} = NaN;
        end

        continue
    end

    %% Read video summarization data
    L_start = cell_BatchVideoSummarization{i,5};
    R_end = cell_BatchVideoSummarization{i,6};
    tFrames = cell_BatchVideoSummarization{i,7};

    L_start_orig = L_start;
    R_end_orig = R_end;
    tFrames_orig = tFrames;

    %% Experimental variables
    [EXP_numpChips,EXP_solvent,EXP_run,EXP_date,EXP_H_channel,EXP_W_channel,EXP_flowrate,EXP_desiredframerate,EXP_collection] = B2KExpParam(vid);
    EXP_actualframerate = vidObj.FrameRate;

    %% Figure Setup
    % REMOVE - set up in function

    %% B2KCoordinateSystem - Image calibration involving channel wall detection and mapping data to appropriate coordinate system
    CS_string = 'CR';
    try
        [I, CS, WR, CR, NR, x_pos_caption, y_pos_caption,myFigures,...
            resizeScale, edgeMethod, R_highest_resolution,T_highest_resolution,...
            H_highest_numpeaks,Hlines_highest_fillgap,Hlines_highest_minlength,...
            Hlines_highest,Hlines_highest_num] = B2KCoordinateSystem(vidObj,L_start_orig, HLcorrectionidx,EXP_W_channel,CS_string,myFigures,flagDebug,flagRunWithDisplay);
    catch ME
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
            ["Video Reader Object" ...
            "Video File Name" ...
            "Pass Count" ...
            "Start Frame" ...
            "Resize Scale" ...
            "Edge Method" ...
            "Hough Rho Resolution" ...
            "Hough Theta Resolution" ...
            "Hough Number of Peaks" ...
            "Houghlines Fill Gap" ...
            "Houghlines Minimum Length" ...
            "Houghlines" ...
            "Houghlines Number Detected"]);

        %% 1) Save table results to .mat
        wallDetTableOutputName = append("table_",grandchildOutputFolderNumberedName,"_wallDet",".mat");
        wallDetTableOutputPath = fullfile(wallDetOutputFolderPath,wallDetTableOutputName);

        % save(fileNamePath);
        save(wallDetTableOutputPath,"table_BatchChannelDetection");

        %% 2) Save used functions to log .txt file

        currScriptNameWithExt = append(currScriptName,'.m');
        [fList,pList] = matlab.codetools.requiredFilesAndProducts(currScriptNameWithExt);

        logFileName = append(grandchildOutputFolderNumberedName,'_log','.txt');
        logFileFullPath = fullfile(wallDetOutputFolderPath,logFileName);
        writelines(grandchildOutputFolderNumberedName,logFileFullPath);

        writelines(currScriptNameWithExt,logFileFullPath,"WriteMode","append");
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
        if i == loopEnd
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
    save(fileNamePath, allvars(tosavefinal).name);

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
        cell_BatchChannelDetection{i,1} = myread(fds.Files{i});
    else
        cell_BatchChannelDetection{i,1} = vidObj;
    end
    cell_BatchChannelDetection{i,2} = cell_BatchChannelDetection{i,1}.Name;
    cell_BatchChannelDetection{i,3} = cell_BatchVideoSummarization{i,5};
    cell_BatchChannelDetection{i,4} = cell_BatchVideoSummarization{i,6};
    cell_BatchChannelDetection{i,5} = resizeScale;
    cell_BatchChannelDetection{i,6} = edgeMethod;
    cell_BatchChannelDetection{i,7} = R_highest_resolution;
    cell_BatchChannelDetection{i,8} = T_highest_resolution;
    cell_BatchChannelDetection{i,9} = H_highest_numpeaks;
    cell_BatchChannelDetection{i,10} = Hlines_highest_fillgap;
    cell_BatchChannelDetection{i,11} = Hlines_highest_minlength;
    cell_BatchChannelDetection{i,12} = Hlines_highest;
    cell_BatchChannelDetection{i,13} = Hlines_highest_num;


    %% Clear figure every loop iteration; close at end of loop
    if i == loopEnd
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

    %% B2KParticleTracking - Detect and track particle

    % Plot variables
    tile_x_pos = -200;
    tile_y_pos = 30;

    %
    try
        [myFigures,position_x_CS,position_y_CS] = B2KParticleTracking03_orig(i,vidObj,myFigures,onePTrackOutputFolderPath,L_start,R_end,tFrames,nFrames,I,CS,x_pos_caption,y_pos_caption,EXP_actualframerate,flagDebug,flagRunWithDisplay);
    catch ME
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

        table_ParticleTracking = cell2table(cell_ParticleTracking,...
            "VariableNames",...
            ["Fill_1" ...
            "Fill_2" ...
            "Fill_3" ...
            "Fill_4" ...
            "Fill_5" ...
            "Fill_6" ...
            "Fill_7" ...
            "Fill_8" ...
            "Fill_9" ...
            "Fill_10" ...
            "Fill_11" ...
            "Fill_12" ...
            "Fill_13"]);

        %% 1) Save table results to .mat
        trackTableOutputName = append("table_",grandchildOutputFolderNumberedName,"_Track",".mat");
        trackTableOutputPath = fullfile(wallDetOutputFolderPath,trackTableOutputName);

        % save(fileNamePath);
        save(trackTableOutputPath,"table_ParticleTracking");

        %% 2) Save used functions to log .txt file

        currScriptNameWithExt = append(currScriptName,'.m');
        [fList,pList] = matlab.codetools.requiredFilesAndProducts(currScriptNameWithExt);

        logFileName = append(grandchildOutputFolderNumberedName,'_log','.txt');
        logFileFullPath = fullfile(wallDetOutputFolderPath,logFileName);
        writelines(grandchildOutputFolderNumberedName,logFileFullPath);

        writelines(currScriptNameWithExt,logFileFullPath,"WriteMode","append");
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
        if i == loopEnd
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

    % -- If no error in B2KParticleTracking, then: --
    % 3) Save workspace variables in .mat
    % 4) Save figure in .tif
    % 5) Store data in cell array

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
    save(fileNamePath, allvars(tosavefinal).name);

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
    cell_ParticleTracking{i,1} = NaN;
    cell_ParticleTracking{i,2} = NaN;
    cell_ParticleTracking{i,3} = NaN;
    cell_ParticleTracking{i,4} = NaN;
    cell_ParticleTracking{i,5} = NaN;
    cell_ParticleTracking{i,6} = NaN;
    cell_ParticleTracking{i,7} = NaN;
    cell_ParticleTracking{i,8} = NaN;
    cell_ParticleTracking{i,9} = NaN;
    cell_ParticleTracking{i,10} = NaN;
    cell_ParticleTracking{i,11} = NaN;
    cell_ParticleTracking{i,12} = NaN;
    cell_ParticleTracking{i,13} = NaN;

    %% Clear figure every loop iteration; close at end of loop
    if i == loopEnd
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

        % flagEndLoopIter = false;

        %% Reset figure counter to initial figure number for loop
        myFigures.figNum = myFigures.initFigNumLoop;
    end

    %%
    %i %report index to keep track of computing time
    i

end

%% Video summarization - log and data table

table_BatchVideoSummarization = cell2table(cell_BatchVideoSummarization,...
    "VariableNames",...
    ["Video Reader Object" ...
    "Video File Name" ...
    "Error Count" ...
    "Pass Count" ...
    "Start Frame" ...
    "End Frame" ...
    "Number of Frames to Analyze" ...
    "Total Number of Frames" ...
    "Flag p-Chip Present" ...
    "Flag p-Chip Stuck" ...
    "Flag no p-Chip" ...
    "Flag Video Corrupt" ...
    "Threshold Median Index" ...
    "Threshold Median Peaks" ...
    "Threshold Median Count" ...
    "Threshold Mean Index" ...
    "Threshold Mean Peaks" ...
    "Threshold Mean Count" ...
    "Threshold Quartiles Index" ...
    "Threshold Quartiles Peaks" ...
    "Threshold Quartiles Count" ...
    "Threshold GESD Index" ...
    "Threshold GESD Peaks" ...
    "Threshold GESD Count" ...
    "TF" ...
    "S1" ...
    "S2" ...
    "ipt" ...
    "residual"...
    "Start StDev" ...
    "End StDev"]);

%% Video summarization - save table results to .mat file

vidSumTableOutputName = append("table_",vidSumFileName,".mat");
vidSumTableOutputPath = fullfile(vidSumOutputFolderPath,vidSumTableOutputName);

% save(fileNamePath);
save(vidSumTableOutputPath,"table_BatchVideoSummarization");

%% Channel detection - log and data table

table_BatchChannelDetection = cell2table(cell_BatchChannelDetection,...
    "VariableNames",...
    ["Video Reader Object" ...
    "Video File Name" ...
    "Pass Count" ...
    "Start Frame" ...
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
save(wallDetTableOutputPath,"table_BatchChannelDetection");


%% Particle tracking - log and data table
table_ParticleTracking = cell2table(cell_ParticleTracking,...
    "VariableNames",...
    ["Fill_1" ...
    "Fill_2" ...
    "Fill_3" ...
    "Fill_4" ...
    "Fill_5" ...
    "Fill_6" ...
    "Fill_7" ...
    "Fill_8" ...
    "Fill_9" ...
    "Fill_10" ...
    "Fill_11" ...
    "Fill_12" ...
    "Fill_13"]);

%% Particle tracking

trackTableOutputName = append("table_",grandchildOutputFolderNumberedName,"_Track",".mat");
trackTableOutputPath = fullfile(wallDetOutputFolderPath,trackTableOutputName);

% save(fileNamePath);
save(trackTableOutputPath,"table_ParticleTracking");

%% Save used functions to log .txt file
currScriptNameWithExt = append(currScriptName,'.m');
[fList,pList] = matlab.codetools.requiredFilesAndProducts(currScriptNameWithExt);

logFileName = append(grandchildOutputFolderNumberedName,'_log','.txt');
logFileFullPath = fullfile(grandchildOutputFolderPath,logFileName);
writelines(grandchildOutputFolderNumberedName,logFileFullPath);

writelines(currScriptNameWithExt,logFileFullPath,"WriteMode","append");
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

%% [FIGURE TEST]
myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);

% Set CreateFcn callback
set(myFigures.fig.(myFigures.handle(myFigures.figNum)), 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 

%
clear B2KMyFigures
% save figure
close(myFigures.fig.(myFigures.handle(myFigures.figNum)))

%% Clear functions
clear functions

%% Timing - whole script, end
tScriptEnd = toc(tScriptStart); % pair 1: toc

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
save(profileOutputFullPath,"profileHandle")
clear profileHandle

% To check Profiler results, 'load' and use 'profview' to load profile handle, Example:
% load("F:\OutputFolder\Batch1pTracking\Batch1pTracking_230816_07\Batch1pTracking_230816_07_profile.mat")
% profview(0,profileHandle)

%% FUNCTIONS - custom reader

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
    home = getenv('USERPROFILE'); if isempty(home), home = char(java.lang.System.getProperty('user.home')); end
    out = fullfile(home,'Downloads');
elseif ismac || isunix
    home = getenv('HOME'); if isempty(home), home = char(java.lang.System.getProperty('user.home')); end
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
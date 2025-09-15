function [settingAddedFrames,flagInternalDrive,flagBatch,flagDebug,flagRunWithDisplay,flagSaveLocDefault,flagTimeFunction,flagTimeSaveFunctionData,flagCalcParameters,flagDataPathOverride,flagDriveOverride,flagSingleVidOverride] = B2KFlags30(addedFrames,drivePriority,processingMode,runMode,displayMode,saveLocation,timeFunction,timeSaveFunctionData,calcParameters,dataPathOverride,driveOverride,singleVidOverride)
%% B2KRunMode - [Function]
%
% Author:  Raymond Yeung
%          Department of Bioengineering
%          University of Riverside, Riverside, CA, 92521, USA
%          ryeun003@ucr.edu
% Release Date: 2025
% B2K Group, Dept. of Bioengineering, Univ. of California, Riverside
% Victor G. J. Rodgers Dept. of Bioengineering, Univ. of California, Riverside
% William H. Grover, Dept. of Bioengineering, Univ. of California, Riverside
% Philip L. Brisk Dept. of Computer Science, Univ. of California, Riverside
%
% Drive priority: internal, external
% Processing modes = single, batch
% Run modes: debug, final
% Display modes: runWithDisplay, runWithoutDisplay
% Save location: defaultVideoLocation, anotherLocation

%% Arguments
arguments
    %---Required---
    addedFrames double
    %---Optional---
    drivePriority {mustBeMember(drivePriority,{'external','internal'})} = 'external'
    processingMode {mustBeMember(processingMode,{'batch','single'})} = 'batch'
    runMode {mustBeMember(runMode,{'debug','final'})} = 'debug'
    displayMode {mustBeMember(displayMode,{'runWithDisplay','runWithoutDisplay'})} = 'runWithDisplay'
    saveLocation {mustBeMember(saveLocation,{'defaultVideoLocation','anotherLocation'})} = 'defaultVideoLocation'
    timeFunction {mustBeMember(timeFunction,{'timeFunction','noTimeFunction'})} = 'timeFunction'
    timeSaveFunctionData {mustBeMember(timeSaveFunctionData,{'timeSaveFunctionData','noTimeSaveFunctionData'})} = 'timeSaveFunctionData'
    calcParameters {mustBeMember(calcParameters,{'calcParameters','noCalcParameters'})} = 'calcParameters'
    dataPathOverride {mustBeText} = ''
    driveOverride {mustBeText} = ''
    singleVidOverride {mustBeText} = ''
    %---Repeating---
end

%% Added Frames
settingAddedFrames = addedFrames;

%% Drive Priority
if matches(drivePriority,'internal')
    flagInternalDrive = true;
else
    flagInternalDrive = false;
end

%% Processing Mode
if matches(processingMode,'batch')
    flagBatch = true;
else
    flagBatch = false;
end

%% Run Mode
if matches(runMode,'debug')
    flagDebug = true;
else
    flagDebug = false;
end

%% Display Mode
if matches(displayMode,'runWithDisplay')
    flagRunWithDisplay = true;
else
    flagRunWithDisplay = false;
end

%% Save Location
if matches(saveLocation,'defaultVideoLocation')
    flagSaveLocDefault = true;
else
    flagSaveLocDefault = false;
end

%% Time Function
if matches(timeFunction,'timeFunction')
    flagTimeFunction = true;
else
    flagTimeFunction = false;
end

%% Time Save Function Data
if matches(timeSaveFunctionData,'timeSaveFunctionData')
    flagTimeSaveFunctionData = true;
else
    flagTimeSaveFunctionData = false;
end

%% Calculation of Parameters
if matches(calcParameters,'calcParameters')
    flagCalcParameters = true;
else
    flagCalcParameters = false;
end

%% Data Path Override
if ~isempty(dataPathOverride)
    flagDataPathOverride = true;
else
    flagDataPathOverride = false;
end

%% Drive Override
if ~isempty(driveOverride)
    flagDriveOverride = true;
else
    flagDriveOverride = false;
end

%% Single Vid Override
if ~isempty(singleVidOverride)
    flagSingleVidOverride = true;
else
    flagSingleVidOverride = false;
end

end
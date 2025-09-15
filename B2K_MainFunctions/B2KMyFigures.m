%% B2KMyFigures.m - [Function] Manage figure numbers and handles using structures containing fields with dynamic expression
%
% Author: Raymond Yeung
% Release date: 2025
% E-mail: ryeun003@ucr.edu
% B2K Group, Dept. of Bioengineering, Univ. of California, Riverside
% Victor G. J. Rodgers Dept. of Bioengineering, Univ. of California, Riverside
% William H. Grover, Dept. of Bioengineering, Univ. of California, Riverside
% Philip L. Brisk Dept. of Computer Science, Univ. of California, Riverside
%
% [STATUS] - incomplete
% [CURRENT FUNC]
% -
% [WORKING ON]
% -
% [BUGS]
% -
% [STILL NEED]
% -
%

function myFigures = B2KMyFigures(myFigures,flagRunWithDisplay)

persistent processNum

if isempty(processNum)
    processNum = 0;
end

processNum = processNum + 1;
myFigures.processNum = processNum;
myFigures.processNumArray = [myFigures.processNumArray;myFigures.processNum];

myFigures.counterNum = myFigures.counterNum + 1;
myFigures.counterNumArray = [myFigures.counterNumArray;myFigures.counterNum];

myFigures.figNum = myFigures.figNum + 1;

% if flagRunWithDisplay == 1
%     myFigures.fig.(B2KFigureDynamicFieldName) = figure(myFigures.figNum,'Visible','on');
% else
%     myFigures.fig.(B2KFigureDynamicFieldName) = figure(myFigures.figNum,'Visible','off');
% end

myFigures.fig.(B2KFigureDynamicFieldName) = figure(myFigures.figNum);
if flagRunWithDisplay == 1
    myFigures.fig.(B2KFigureDynamicFieldName).Visible = 'on';
else
    myFigures.fig.(B2KFigureDynamicFieldName).Visible = 'off';
end


myFigures.handle(myFigures.figNum,1) = B2KFigureDynamicFieldName;
myFigures.figNumArray = [myFigures.figNumArray;myFigures.figNum];

    function dynFieldName = B2KFigureDynamicFieldName
        formatSpec = "f%d";
        dynFieldName = sprintf(formatSpec,myFigures.figNum);
    end

% Set monitor which the figure displays on
monitorInfo = get(0,"MonitorPositions");
% Find highest resolution monitor and use to display figure
[monitorResScaled,~] = max(monitorInfo(:,3)); %finds only first occurence
maxIdx = find(monitorInfo==monitorResScaled); %used to find all occurrences
szMonitorInfo = size(monitorInfo);
[monitorNum,~] = ind2sub(szMonitorInfo,maxIdx);

if size(monitorNum,1) > 1 % If more than one monitor with max res, then use former
    monitorChoice = max(monitorNum);
else
    monitorChoice = monitorNum; % If just one monitor with max res, then use that one
end

myFigures.fig.(myFigures.handle(myFigures.figNum)).Position = monitorInfo(monitorChoice,:);
myFigures.fig.(myFigures.handle(myFigures.figNum)).WindowState ='maximized';

end
function [L_start,R_end,tFrames,nFrames,...
    flag_pChip_present,flag_pChip_stuck,flag_no_pChip,flag_video_corrupt,...
    myFigures,varargout] = B2KVideoSummary01(vidObj,myFigures,flagDebug,flagRunWithDisplay,addedFrames)
%% B2KVideoSummary - [Function] Summarizes range of video frames that consist of objects
%
%
% Author:  Raymond K. Yeung
%          Department of Bioengineering
%          University of Riverside, Riverside, CA, 92521, USA
%          ryeun003@ucr.edu
% Release Date: 2023
% B2K Group, Dept. of Bioengineering, Univ. of California, Riverside
% Victor G. J. Rodgers Dept. of Bioengineering, Univ. of California, Riverside
% William H. Grover, Dept. of Bioengineering, Univ. of California, Riverside
% Philip L. Brisk Dept. of Computer Science, Univ. of California, Riverside
% [Notes]
% 

%% Inputs
% Required
%   vidObj
%   flagDebug
%   flagRunWithDisplay

% Optional
%   addedFrames

% Name-Value

%% Arguments
arguments
    vidObj (1,1) VideoReader
    myFigures (1,1) struct
    flagDebug (1,1) logical
    flagRunWithDisplay (1,1) logical
    addedFrames (1,1) double {mustBeNumeric}
end

%% Final Mode

% Error flags
flag_no_pChip = NaN;
flag_pChip_present = NaN;
flag_pChip_stuck = NaN;
flag_video_corrupt = NaN;

% Preallocate arrays
L_start = NaN;
R_end = NaN;
tFrames = NaN;

nFrames = vidObj.NumFrames;
stdImage = NaN(nFrames,1);

% Handle corrupt videos - return empty figure and NaN values
if nFrames == 0
    flag_video_corrupt = 1;

    myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);
    thresh_median_idx = NaN;
    thresh_median_pks = NaN;
    thresh_median_count = NaN;
    thresh_mean_idx = NaN;
    thresh_mean_pks = NaN;
    thresh_mean_count = NaN;
    thresh_quartiles_idx = NaN;
    thresh_quartiles_pks = NaN;
    thresh_quartiles_count = NaN;
    thresh_gesd_idx = NaN;
    thresh_gesd_pks = NaN;
    thresh_gesd_count = NaN;
    TF = NaN;
    S1 = NaN;
    S2 = NaN;
    ipt = NaN;
    residual = NaN;
    L_startBlockStd = NaN;
    R_endBlockStd = NaN;
    return
end

%% Number of tiles
numTiles = 9;


%% Parameter controlling frames to add before and after evaluated endpoints
% addedFrames = 25;

%% Debug Mode

%% Run with Display

%% Run Without Display


%% Preallocate arrays


%% Evaluate standard deviation of image for every frame

for i = 1:nFrames
    fr = read(vidObj,i);
    fr_bw = rgb2gray(fr);
    stdImage(i) = std2(fr_bw);
end

stdImage_smoothed = smoothdata(stdImage,"movmean",25); %still fixed #479

stdImage_smoothed_d = gradient(stdImage_smoothed);
frame = 1:1:length(stdImage_smoothed);

% Smoothed SD of image
[outlier_median_sm,outlier_median_L_sm,outlier_median_U_sm,outlier_median_C_sm] = isoutlier(stdImage_smoothed,"median");
[outlier_mean_sm,outlier_mean_L_sm,outlier_mean_U_sm,outlier_mean_C_sm] = isoutlier(stdImage_smoothed,"mean");
[outlier_quartiles_sm,outlier_quartiles_L_sm,outlier_quartiles_U_sm,outlier_quartiles_C_sm] = isoutlier(stdImage_smoothed,"quartiles");
[outlier_gesd_sm,outlier_gesd_L_sm,outlier_gesd_U_sm,outlier_gesd_C_sm] = isoutlier(stdImage_smoothed,"gesd");

% Smoothed 1st Derivative of SD of image
[outlier_median_sm_d,outlier_median_L_sm_d,outlier_median_U_sm_d,outlier_median_C_sm_d] = isoutlier(stdImage_smoothed_d,"median");
[outlier_mean_sm_d,outlier_mean_L_sm_d,outlier_mean_U_sm_d,outlier_mean_C_sm_d] = isoutlier(stdImage_smoothed_d,"mean");
[outlier_quartiles_sm_d,outlier_quartiles_L_sm_d,outlier_quartiles_U_sm_d,outlier_quartiles_C_sm_d] = isoutlier(stdImage_smoothed_d,"quartiles");
[outlier_gesd_sm_d,outlier_gesd_L_sm_d,outlier_gesd_U_sm_d,outlier_gesd_C_sm_d] = isoutlier(stdImage_smoothed_d,"gesd",'ThresholdFactor',0.05); %ThresholdFactor default 0.05

%% Create figure
myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);

%% Statistics to determine appropriate threshold(s)

%MinPeakHeight - Minimum peak height
%MinPeakProminence - Minimum peak prominence
%Threshold - Minimum height difference
%MinPeakDistance - Minimum peak separation
%MinPeakWidth - Minimum peak width
%MaxPeakWidth - Maximum peak width

%WidthReference - Reference height for width measurements

% Conditions for finding peaks
valMinPeakProminence = 0.0065;
valMinPeakWidth = 0;

[thresh_median_pos_pks,thresh_median_pos_locs,thresh_median_pos_w,thresh_median_pos_p] = findpeaks(stdImage_smoothed_d,'MinPeakHeight',outlier_median_U_sm_d,'MinPeakProminence',valMinPeakProminence,'MinPeakWidth',valMinPeakWidth);
[thresh_median_neg_pks,thresh_median_neg_locs,thresh_median_neg_w,thresh_median_neg_p] = findpeaks(-stdImage_smoothed_d,'MinPeakHeight',outlier_median_U_sm_d,'MinPeakProminence',valMinPeakProminence,'MinPeakWidth',valMinPeakWidth);
thresh_median_idx = [thresh_median_pos_locs; thresh_median_neg_locs];
thresh_median_pks = [thresh_median_pos_pks; thresh_median_neg_pks];
thresh_median_count = length(thresh_median_idx);

[thresh_mean_pos_pks,thresh_mean_pos_locs,thresh_mean_pos_w,thresh_mean_pos_p] = findpeaks(stdImage_smoothed_d,'MinPeakHeight',outlier_mean_U_sm_d,'MinPeakProminence',valMinPeakProminence,'MinPeakWidth',valMinPeakWidth);
[thresh_mean_neg_pks,thresh_mean_neg_locs,thresh_mean_neg_w,thresh_mean_neg_p] = findpeaks(-stdImage_smoothed_d,'MinPeakHeight',outlier_mean_U_sm_d,'MinPeakProminence',valMinPeakProminence,'MinPeakWidth',valMinPeakWidth);
thresh_mean_idx = [thresh_mean_pos_locs; thresh_mean_neg_locs];
thresh_mean_pks = [thresh_mean_pos_pks; thresh_mean_neg_pks];
thresh_mean_count = length(thresh_mean_idx);

[thresh_quartiles_pos_pks,thresh_quartiles_pos_locs,thresh_quartiles_pos_w,thresh_quartiles_pos_p] = findpeaks(stdImage_smoothed_d,'MinPeakHeight',outlier_quartiles_U_sm_d,'MinPeakProminence',valMinPeakProminence,'MinPeakWidth',valMinPeakWidth);
[thresh_quartiles_neg_pks,thresh_quartiles_neg_locs,thresh_quartiles_neg_w,thresh_quartiles_neg_p] = findpeaks(-stdImage_smoothed_d,'MinPeakHeight',outlier_quartiles_U_sm_d,'MinPeakProminence',valMinPeakProminence,'MinPeakWidth',valMinPeakWidth);
thresh_quartiles_idx = [thresh_quartiles_pos_locs; thresh_quartiles_neg_locs];
thresh_quartiles_pks = [thresh_quartiles_pos_pks; thresh_quartiles_neg_pks];
thresh_quartiles_count = length(thresh_quartiles_idx);

[thresh_gesd_pos_pks,thresh_gesd_pos_locs,thresh_gesd_pos_w,thresh_gesd_pos_p] = findpeaks(stdImage_smoothed_d,'MinPeakHeight',outlier_gesd_U_sm_d,'MinPeakProminence',valMinPeakProminence,'MinPeakWidth',valMinPeakWidth);
[thresh_gesd_neg_pks,thresh_gesd_neg_locs,thresh_gesd_neg_w,thresh_gesd_neg_p] = findpeaks(-stdImage_smoothed_d,'MinPeakHeight',outlier_gesd_U_sm_d,'MinPeakProminence',valMinPeakProminence,'MinPeakWidth',valMinPeakWidth);
thresh_gesd_idx = [thresh_gesd_pos_locs; thresh_gesd_neg_locs];
thresh_gesd_pks = [thresh_gesd_pos_pks; thresh_gesd_neg_pks];
thresh_gesd_count = length(thresh_gesd_idx);

%% Find change points

% Maximum number of change points
MaxNumChanges = 2;

[ipt,residual] = findchangepts(stdImage_smoothed,'MaxNumChanges',MaxNumChanges,'Statistic','mean');

%% Find the start frame and end frame based on different scenarios

whilestart = 1;
while whilestart == 1

    if isempty(thresh_gesd_pos_locs) && isempty(thresh_gesd_neg_locs) % no p-Chip
        L_start_orig = NaN;
        L_start = NaN;
        R_end_orig = NaN;
        R_end = NaN;
        flag_no_pChip = 1;
        break
    elseif isempty(thresh_gesd_pos_locs) && ~isempty(thresh_gesd_neg_locs) % p-Chip present and left
        L_start_orig = 1;
        L_start = 1;
        R_end_orig = thresh_gesd_neg_locs(end);
        R_end = R_end_orig + addedFrames;
        if R_end > nFrames
            R_end = nFrames;
        end
        flag_pChip_present = 1;
        % need another condition for: p-Chip present and stuck
    elseif ~isempty(thresh_gesd_pos_locs) && isempty(thresh_gesd_neg_locs) % p-Chip stuck
        L_start_orig = thresh_gesd_pos_locs(1);
        L_start = L_start_orig - addedFrames;
        if L_start < 1
            L_start = 1;
        end
        R_end_orig = nFrames;
        R_end = R_end_orig;
        flag_pChip_stuck = 1;
    else % normal
        L_start_orig = thresh_gesd_pos_locs(1);
        L_start = L_start_orig - addedFrames;
        if L_start < 1
            L_start = 1;
        end
        R_end_orig = thresh_gesd_neg_locs(end);
        R_end = R_end_orig + addedFrames;
        if R_end > nFrames
            R_end = nFrames;
        end
        flag_pChip_present = 1;
    end

    % total number of frames to analyze
    tFrames = R_end - L_start + 1;
    whilestart = 0;

end

if flag_no_pChip ~= 1
    I_start = read(vidObj,L_start);
    I_start = rgb2gray(I_start);

    I_mid = read(vidObj,floor((L_start+R_end)/2));
    I_mid = rgb2gray(I_mid);

    I_end = read(vidObj,R_end);
    I_end = rgb2gray(I_end);
end

%%
blockFrames = 50;
if ~isnan(L_start_orig)

    L_startLeftIdx = L_start_orig - blockFrames;
    if L_startLeftIdx < 1
        L_startLeftIdx = 1;
    end

    L_startRightIdx = L_start + blockFrames;
    if L_startRightIdx > nFrames
        L_startRightIdx = nFrames;
    end

    L_startBlock = stdImage_smoothed(L_startLeftIdx:L_startRightIdx);
    L_startBlockStd = std(L_startBlock);
else
    L_startBlockStd = NaN;
end

if ~isnan(R_end_orig)

    R_endLeftIdx = R_end_orig-blockFrames;
    if R_endLeftIdx < 1
        R_endLeftIdx = 1;
    end

    R_endRightIdx = R_end_orig+blockFrames;
    if R_endRightIdx > nFrames
        R_endRightIdx = nFrames;
    end

    R_endBlock = stdImage_smoothed(R_endLeftIdx:R_endRightIdx);
    R_endBlockStd = std(R_endBlock);
else
    R_endBlockStd = NaN;
end

%%

% Set monitor which the figure displays on
monitorInfo = get(0,"MonitorPositions");
% Find highest resolution monitor and use to display figure
[monitorResScaled,~] = max(monitorInfo(:,3)); %finds only first occurence
maxIdx = find(monitorInfo==monitorResScaled); %used to find all occurrences
szMonitorInfo = size(monitorInfo);
[monitorNum,~] = ind2sub(szMonitorInfo,maxIdx);


if size(monitorNum,1) > 1 % If more than one monitor, then use latter
    monitorChoice = max(monitorNum);
else
    monitorChoice = monitorNum; % If more than one monitor, then use former
end

myFigures.fig.(myFigures.handle(myFigures.figNum)).Position = monitorInfo(monitorChoice,:);
myFigures.fig.(myFigures.handle(myFigures.figNum)).WindowState ='maximized';


tiledlayout(numTiles,1,"TileSpacing","tight","Padding","tight");

% nexttile
% plot(frame,stdImage_smoothed)
% xlabel('Frame # [-]')
% ylabel({'Smoothed SD of image [-]'})
% yline([outlier_median_U_sm],'--',{'outlier median'},'Color',"#0072BD")
% yline([outlier_mean_U_sm],'--',{'outlier mean'},'Color',"#D95319") %bad
% yline([outlier_quartiles_U_sm],'--b',{'outlier quartiles'},'Color',"#EDB120")
% yline([outlier_gesd_U_sm],'--',{'outlier gesd'},'Color',"#7E2F8E") %bad
% yline([outlier_median_L_sm],'--',{'outlier median'},'Color',"#0072BD") %good
% yline([outlier_mean_L_sm],'--',{'outlier mean'},'Color',"#D95319") %bad
% yline([outlier_quartiles_L_sm],'--',{'outlier quartiles'},'Color',"#EDB120") %better
% yline([outlier_gesd_L_sm],'--',{'outlier gesd'},'Color',"#7E2F8E") %bad
% yline([mean(stdImage_smoothed)-std2(stdImage_smoothed),mean(stdImage_smoothed)+std2(stdImage_smoothed)],'--',{'1\sigma','1\sigma'},'Color',"#77AC30")
% yline([mean(stdImage_smoothed)-2*std2(stdImage_smoothed),mean(stdImage_smoothed)+2*std2(stdImage_smoothed)],'--',{'2\sigma','2\sigma'},'Color',"#4DBEEE")
% yline([mean(stdImage_smoothed)-3*std2(stdImage_smoothed),mean(stdImage_smoothed)+3*std2(stdImage_smoothed)],'--',{'3\sigma','3\sigma'},'Color',"#A2142F")

%%
nexttile

[TF,S1,S2] = ischange(stdImage_smoothed,'MaxNumChanges',MaxNumChanges);

plot(frame,stdImage_smoothed)

hold on
stairs(S1,'Color','r','LineWidth',1)

xlabel('Frame # [-]')
formatSpec = 'Total res. error: %.4f';
str_start = sprintf(formatSpec,residual);
ylabel({'Smoothed SD of image [-]';str_start})
yline([outlier_median_U_sm],'--',{'outlier median'},'Color',"#0072BD")
yline([outlier_mean_U_sm],'--',{'outlier mean'},'Color',"#D95319") %bad
yline([outlier_quartiles_U_sm],'--b',{'outlier quartiles'},'Color',"#EDB120")
yline([outlier_gesd_U_sm],'--',{'outlier gesd'},'Color',"#7E2F8E") %bad
yline([outlier_median_L_sm],'--',{'outlier median'},'Color',"#0072BD") %good
yline([outlier_mean_L_sm],'--',{'outlier mean'},'Color',"#D95319") %bad
yline([outlier_quartiles_L_sm],'--',{'outlier quartiles'},'Color',"#EDB120") %better
yline([outlier_gesd_L_sm],'--',{'outlier gesd'},'Color',"#7E2F8E") %bad
yline([mean(stdImage_smoothed)-std2(stdImage_smoothed),mean(stdImage_smoothed)+std2(stdImage_smoothed)],'--',{'1\sigma','1\sigma'},'Color',"#77AC30")
yline([mean(stdImage_smoothed)-2*std2(stdImage_smoothed),mean(stdImage_smoothed)+2*std2(stdImage_smoothed)],'--',{'2\sigma','2\sigma'},'Color',"#4DBEEE")
yline([mean(stdImage_smoothed)-3*std2(stdImage_smoothed),mean(stdImage_smoothed)+3*std2(stdImage_smoothed)],'--',{'3\sigma','3\sigma'},'Color',"#A2142F")


hold off

%% Area evaluation
nexttile

stdImage_smoothed_n = stdImage_smoothed - min(S1);

outlier_median_U_sm_norm = outlier_median_U_sm - min(S1);
outlier_mean_U_sm_norm = outlier_mean_U_sm - min(S1);
outlier_quartiles_U_sm_norm = outlier_quartiles_U_sm - min(S1);
outlier_gesd_U_sm_norm = outlier_gesd_U_sm - min(S1);

plot(frame,stdImage_smoothed_n)
hold on

xlabel('Frame # [-]')
ylabel({'Smoothed SD of image';'with normalization [-]'});
yline([outlier_median_U_sm_norm],'--',{'outlier median'},'Color',"#0072BD")
yline([outlier_mean_U_sm_norm],'--',{'outlier mean'},'Color',"#D95319") %bad
yline([outlier_quartiles_U_sm_norm],'--b',{'outlier quartiles'},'Color',"#EDB120")
yline([outlier_gesd_U_sm_norm],'--',{'outlier gesd'},'Color',"#7E2F8E") %bad
yline([mean(stdImage_smoothed)+std2(stdImage_smoothed)-min(S1)],'--',{'1\sigma'},'Color',"#77AC30")
yline([mean(stdImage_smoothed)+2*std2(stdImage_smoothed)-min(S1)],'--',{'2\sigma'},'Color',"#4DBEEE")
yline([mean(stdImage_smoothed)+3*std2(stdImage_smoothed)-min(S1)],'--',{'3\sigma'},'Color',"#A2142F")

% Threshold for height
% SDpkThreshold = 0.075;
SDpkThreshold = 0.050;

[SDpks,SDlocs,SDw,SDp] = findpeaks(stdImage_smoothed_n,'MinPeakHeight',SDpkThreshold);



peakInd = find(stdImage_smoothed_n > SDpkThreshold);
if ~isempty(peakInd)

    grouped = mat2cell(peakInd.', 1, diff( [0, find(diff(peakInd.') ~= 1), length(peakInd.')] )) ;

    plot(SDlocs,SDpks,'o')

    groupedLength = length(grouped);
    area_peak = zeros(1,groupedLength);

    % Set area threshold condition
    area_thresh = 25;

    for k = 1:groupedLength
        groupedIdx = grouped{k};
        area_peak(k) = trapz(groupedIdx,stdImage_smoothed_n(groupedIdx));

        plot(groupedIdx,stdImage_smoothed_n(groupedIdx),'r')
        text(frame(groupedIdx(end)),stdImage_smoothed_n(groupedIdx(end)), compose('\\leftarrow Area = %.2f', area_peak(k)), 'Horiz','left', 'Vert','middle', 'Rotation',30, 'FontSize',10)

        % Color
        if area_peak(k) > area_thresh
            plot(groupedIdx,stdImage_smoothed_n(groupedIdx),'g')
        end
    end

    % Filter: remove peak(s) with area < 25
    groupedAreaFiltered = {};
    groupedAreaCounter = 0;
    for j = 1:groupedLength
        if area_peak(j) >= area_thresh
            groupedAreaCounter = groupedAreaCounter + 1;
            groupedAreaFiltered(:,groupedAreaCounter) = grouped(1,j);
        end
    end

    if length(groupedAreaFiltered) == 1
        L_start_area = groupedAreaFiltered{1,1}(1,1);
        R_end_area = groupedAreaFiltered{1,1}(1,end);
        xline(L_start_area,'--m','pos')
        xline(R_end_area,'--r','neg')
    else
        L_start_area = groupedAreaFiltered{1,1}(1,1);
        R_end_area = groupedAreaFiltered{1,end}(1,end);
        xline(groupedAreaFiltered{1,1}(1,1),'--m','pos')
        xline(groupedAreaFiltered{1,end}(1,end),'--r','neg')
    end

    if stdImage_smoothed_d(L_start_area) < 0.004 %0.01, 0.005, 0.003 (#470 error),
        xline(L_start,'--b','start')
    end

    if -stdImage_smoothed_d(R_end_area) < 0.004 %0.01, 0.005, 0.003 (#470 error),
        xline(R_end,'--g','end')
    end

end

hold off


%%
nexttile
plot(frame,stdImage_smoothed_d)
xlabel('Frame # [-]')
ylabel({'Derivative of'; 'smoothed SD of image [-]'})
yline([outlier_median_U_sm_d],'--',{'outlier median'},'Color',"#0072BD")
yline([outlier_mean_U_sm_d],'--',{'outlier mean'},'Color',"#D95319") %bad
yline([outlier_quartiles_U_sm_d],'--b',{'outlier quartiles'},'Color',"#EDB120")
yline([outlier_gesd_U_sm_d],'--',{'outlier gesd'},'Color',"#7E2F8E") %bad
yline([outlier_median_L_sm_d],'--',{'outlier median'},'Color',"#0072BD") %good
yline([outlier_mean_L_sm_d],'--',{'outlier mean'},'Color',"#D95319") %bad
yline([outlier_quartiles_L_sm_d],'--',{'outlier quartiles'},'Color',"#EDB120") %better
yline([outlier_gesd_L_sm_d],'--',{'outlier gesd'},'Color',"#7E2F8E") %bad
yline([-std2(stdImage_smoothed_d),std2(stdImage_smoothed_d)],'--',{'1\sigma','1\sigma'},'Color',"#77AC30")
yline([-2*std2(stdImage_smoothed_d),2*std2(stdImage_smoothed_d)],'--',{'2\sigma','2\sigma'},'Color',"#4DBEEE")
yline([-3*std2(stdImage_smoothed_d),3*std2(stdImage_smoothed_d)],'--',{'3\sigma','3\sigma'},'Color',"#A2142F")

if flag_no_pChip == 1
    nexttile
    imshow(read(vidObj,1))
    ylabel('No-Chip')

    nexttile
    imshow(read(vidObj,floor((1+nFrames)/2)))
    ylabel('No p-Chip')

    nexttile
    imshow(read(vidObj,nFrames))
    ylabel('No p-Chip')
elseif flag_pChip_present == 1 && isnan(flag_pChip_stuck)
    xline(L_start,'--m','pos')
    xline(R_end,'--r','neg')

    nexttile
    imshow(read(vidObj,1))
    ylabel('p-Chip present and left')

    nexttile
    imshow(read(vidObj,floor((1+nFrames)/2)))
    ylabel('p-Chip present and left')

    nexttile
    imshow(read(vidObj,nFrames))
    ylabel('p-Chip present and left')
elseif flag_pChip_present == 1 && flag_pChip_stuck == 1
    xline(L_start,'--m','pos')
    xline(R_end,'--r','neg')

    nexttile
    imshow(read(vidObj,1))
    ylabel('p-Chip present and stuck')

    nexttile
    imshow(read(vidObj,floor((1+nFrames)/2)))
    ylabel('p-Chip present and stuck')

    nexttile
    imshow(read(vidObj,nFrames))
    ylabel('p-Chip present and stuck')
elseif isnan(flag_pChip_present) && flag_pChip_stuck == 1
    xline(L_start,'--m','pos')
    xline(R_end,'--r','neg')

    nexttile
    imshow(read(vidObj,1))
    ylabel('p-Chip stuck')

    nexttile
    imshow(read(vidObj,floor((1+nFrames)/2)))
    ylabel('p-Chip stuck')

    nexttile
    imshow(read(vidObj,nFrames))
    ylabel('p-Chip stuck')
else
    xline(L_start,'--m','pos')
    xline(R_end,'--r','neg')

    nexttile
    imshow(I_start)
    formatSpec = 'Start frame: %.f';
    str_start = sprintf(formatSpec,L_start_orig);
    ylabel(str_start);

    nexttile
    imshow(I_mid)
    formatSpec = 'Middle frame: %.f';
    str_start = sprintf(formatSpec,floor((L_start_orig+R_end_orig)/2));
    ylabel(str_start);

    nexttile
    imshow(I_end)
    formatSpec = 'End frame: %.f';
    str_end = sprintf(formatSpec,R_end_orig);
    ylabel(str_end);
end

%% Regionprops

if isnan(L_start_orig)
    num_fr_start = 1;
else
    num_fr_start = L_start_orig;
end

if isnan(R_end_orig)
    num_fr_end = nFrames;
else
    num_fr_end = R_end_orig;
end

nexttile
fr_start = read(vidObj,num_fr_start);
fr_start_gray = rgb2gray(fr_start);
fr_start_bw = ~imbinarize(fr_start_gray,0.32);
% fr_start_bw_areaopen = bwareaopen(fr_start_bw,50);
imshow(fr_start_bw)
stats_start = regionprops(fr_start_bw,"Area","BoundingBox","Centroid","Circularity");

% stats_start_pChip = struct([]);
stats_start_pChip_counter = 0;

for k = 1:length(stats_start)
    if stats_start(k).Area >= 50 && stats_start(k).Circularity >= 0.5
        stats_start_pChip_counter = stats_start_pChip_counter + 1;
        stats_start_pChip(stats_start_pChip_counter) = stats_start(k);
    end
end

% fr_start_bw2 = bwpropfilt(fr_start_bw,'Area',[50 500]);
% imshow(fr_start_bw2)



nexttile
fr_mid = read(vidObj,floor((num_fr_start+num_fr_end)/2));
fr_mid_gray = rgb2gray(fr_mid);
fr_mid_bw = ~imbinarize(fr_mid_gray,0.32);
% fr_end_bw_areaopen = bwareaopen(fr_end_bw,50);
imshow(fr_mid_bw)
stats_end = regionprops(fr_mid_bw,"Area","BoundingBox","Centroid","Circularity");

% stats_end_pChip = struct([]);
stats_end_pChip_counter = 0;

for k = 1:length(stats_end)
    if stats_end(k).Area >= 50 && stats_end(k).Circularity >= 0.5
        stats_end_pChip_counter = stats_end_pChip_counter + 1;
        stats_end_pChip(stats_end_pChip_counter) = stats_end(k);
    end
end



nexttile
fr_end = read(vidObj,nFrames);
fr_end_gray = rgb2gray(fr_end);
fr_end_bw = ~imbinarize(fr_end_gray,0.32);
% fr_end_bw_areaopen = bwareaopen(fr_end_bw,50);
imshow(fr_end_bw)
stats_end = regionprops(fr_end_bw,"Area","BoundingBox","Centroid","Circularity");

% stats_end_pChip = struct([]);
stats_end_pChip_counter = 0;

for k = 1:length(stats_end)
    if stats_end(k).Area >= 50 && stats_end(k).Circularity >= 0.5
        stats_end_pChip_counter = stats_end_pChip_counter + 1;
        stats_end_pChip(stats_end_pChip_counter) = stats_end(k);
    end
end

%% Variable output arguments
varargout{1} = thresh_median_idx;
varargout{2} = thresh_median_pks;
varargout{3} = thresh_median_count;
varargout{4} = thresh_mean_idx;
varargout{5} = thresh_mean_pks;
varargout{6} = thresh_mean_count;
varargout{7} = thresh_quartiles_idx;
varargout{8} = thresh_quartiles_pks;
varargout{9} = thresh_quartiles_count;
varargout{10} = thresh_gesd_idx;
varargout{11} = thresh_gesd_pks;
varargout{12} = thresh_gesd_count;
varargout{13} = TF;
varargout{14} = S1;
varargout{15} = S2;
varargout{16} = ipt;
varargout{17} = residual;
varargout{18} = L_startBlockStd;
varargout{19} = R_endBlockStd;

%%
clear B2KMyFigures

end


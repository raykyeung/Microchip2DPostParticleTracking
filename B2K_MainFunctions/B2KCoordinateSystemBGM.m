function [I_average, CS, WR, CR, NR, x_pos_caption, y_pos_caption,myFigures,varargout] = B2KCoordinateSystemBGM(vidObj,L_start_orig, HLcorrectionidx,EXP_W_channel,CS_string,myFigures,flagDebug,flagRunWithDisplay)
%% B2KCoordinateSystemBGM.m - [Function]
%
% Author: Raymond Yeung
% Release date: 2025
% E-mail: ryeun003@ucr.edu
% B2K Group, Dept. of Bioengineering, Univ. of California, Riverside
% Victor G. J. Rodgers Dept. of Bioengineering, Univ. of California, Riverside
% William H. Grover, Dept. of Bioengineering, Univ. of California, Riverside
% Philip L. Brisk Dept. of Computer Science, Univ. of California, Riverside

%% Plot variables
% tile_x_pos = -200;
% tile_y_pos = 46;
tile_x_pos = -250;
tile_y_pos = 46;

tile_x_pos_norm = -0.05;
% tile_y_pos_norm = ;

fontSizeYLabel = 10;
txtFontSize = 14;

txtFontSizePixel = 20;

%% Mapping data to coordinate systems
while 1 %check for appropriate houghline conditions
    
    %% 1) I
    I = read(vidObj,L_start_orig - HLcorrectionidx);
    I = rgb2gray(I);

    %% 2) Edge detection, I
    edgeMethod = "sobel";
    edge_I = edge(I,edgeMethod);

    %% 3) I_resize
    
    resizeScale = 1.5;
    I_resize = imresize(I,resizeScale);

    %% 4) Edge detection, I_resize
    
    edge_I_resize = edge(I_resize,edgeMethod);

    %% 5) I_average

    nFrames = vidObj.NumFrames;

    for f = 1:nFrames
        if f == 1
            sumFrame = double(rgb2gray(read(vidObj,f)));
        else
            sumFrame = sumFrame + double(rgb2gray(read(vidObj,f)));
        end
    end
    averageFrame = sumFrame / nFrames;
    I_average = uint8(averageFrame);

    %% 6) Edge detection, I_average

    edge_I_average = edge(I_average,edgeMethod);


    %% 7) Houghline, custom - find indices of highest voted in Hough accumulator array
    
    % Hough Parameters
    R_resolution = 0.1;
    T_resolution = 0.1;

    H_numpeaks = 20;
    Hlines_fillgap = 200;
    Hlines_minlength = 300;
%     H_numpeaks = 50;
%     Hlines_fillgap = 1000;
%     Hlines_minlength = 50;

    % Limit theta ranges to plus/minus angle tolerance from horizontal (zero degree angle) line
    Hlines_angletol_pre = 2;

    T_start_angle_negzero_pre = -90; %accounts for zero angle
    T_end_angle_negzero_pre = -90 + Hlines_angletol_pre;
    [H1_pre,T1_pre,R_pre] = hough(edge_I_average,'RhoResolution',R_resolution,'Theta',T_start_angle_negzero_pre:T_resolution:T_end_angle_negzero_pre);

    T_start_angle_pos_pre = 90 - Hlines_angletol_pre; %already accounted for zero degree angle beforehand
    T_end_angle_pos_pre = 90 - T_resolution;
    [H2_pre,T2_pre,R2_pre] = hough(edge_I_average,'RhoResolution',R_resolution,'Theta',T_start_angle_pos_pre:T_resolution:T_end_angle_pos_pre);
    
    H_pre = [H1_pre H2_pre];
    T_pre = [T1_pre T2_pre];

    [sortedH_pre,sortedHidx_pre] = sort(H_pre(:),"descend");
    sortedP_pre = [sortedHidx_pre,sortedH_pre];
    topHvalue_pre = H_numpeaks;

    highestP_pre = sortedP_pre(1:topHvalue_pre,1:2);

    Rsize_pre = size(H_pre,1);
    Tsize_pre = size(H_pre,2);

    highestHlinidx_pre = highestP_pre(:,2).';

    [highestR_pre,highestT_pre] = ind2sub([Rsize_pre Tsize_pre],highestHlinidx_pre);

    P_pre = zeros(topHvalue_pre,2);
    P_pre(:,1) = sortedP_pre(1:topHvalue_pre,1);
    P_pre(:,2) = highestT_pre.';

    P_filter_pre = P_pre(1:topHvalue_pre,:);

    % FillGap and MinLength Filtering
    Hlines_pre = B2Khoughlines(edge_I_average,T_pre,R_pre,P_filter_pre,'FillGap',Hlines_fillgap,'MinLength',Hlines_minlength);
    Hlines_num_pre = length(Hlines_pre);

    %% 8) Rotate and crop image if a houghline angle offset is present
    
    if Hlines_pre(1).point1 >= 0
        detectedAngle = Hlines_pre(1).theta + 90;
    else
        detectedAngle = Hlines_pre(1).theta - 90;
    end

    if detectedAngle == 0
        I_rot = I_average;
        I_rotcrop = I_average;
    else
        I_rot = imrotate(I_average,detectedAngle,'nearest','loose');
        [CI,~] = B2KRotateandCrop(I_average,detectedAngle);
        I_rotcrop = CI;
    end

    %% 9) Edge detection, I_rotcrop

    edge_I_rotcrop = edge(I_rotcrop,edgeMethod);

    %% 10) Houghlines, I_rotcrop

    % Limit theta ranges to plus/minus angle tolerance from horizontal (zero degree angle) line
    Hlines_highest_angletol = 2;

    T_start_angle_negzero = -90; %accounts for zero angle
    T_end_angle_negzero = -90 + Hlines_highest_angletol;
    [H1,T1,R] = hough(edge_I_average,'RhoResolution',R_resolution,'Theta',T_start_angle_negzero:T_resolution:T_end_angle_negzero);

    T_start_angle_pos = 90 - Hlines_highest_angletol; %already accounted for zero degree angle beforehand
    T_end_angle_pos = 90 - T_resolution;
    [H2,T2,R2] = hough(edge_I_average,'RhoResolution',R_resolution,'Theta',T_start_angle_pos:T_resolution:T_end_angle_pos);
    
    H = [H1 H2];
    T = [T1 T2];

    [sortedH,sortedHidx] = sort(H(:),"descend");
    sortedP = [sortedHidx,sortedH];
    topHvalue = H_numpeaks;

    highestP = sortedP(1:topHvalue,1:2);

    Rsize = size(H,1);
    Tsize = size(H,2);

    highestHlinidx = highestP(:,2).';

    [highestR,highestT] = ind2sub([Rsize Tsize],highestHlinidx);

    P = zeros(topHvalue,2);
    P(:,1) = sortedP(1:topHvalue,1);
    P(:,2) = highestT.';

    P_filter = P(1:topHvalue,:);

    % FillGap and MinLength Filtering
    Hlines_highest = B2Khoughlines(edge_I_average,T,R,P_filter,'FillGap',Hlines_fillgap,'MinLength',Hlines_minlength);
    Hlines_highest_num = length(Hlines_highest);

    %% 11) Grouping Houghlines

    % Lines
    Hlines_highest_ypos = NaN(Hlines_highest_num,3);

    for Hlines_highest_ypos_idx = 1:Hlines_highest_num
        Hlines_highest_ypos(Hlines_highest_ypos_idx,1) = Hlines_highest_ypos_idx;
        Hlines_highest_ypos(Hlines_highest_ypos_idx,2) = [Hlines_highest(Hlines_highest_ypos_idx).point1(1,2)];
        Hlines_highest_ypos(Hlines_highest_ypos_idx,3) = [Hlines_highest(Hlines_highest_ypos_idx).votes];
    end

    % Grouping lines based on line divide located at midpoint between lines with the largest difference in y-distance
    [distSortedHlines,distSortedHlines_idx] = sort(Hlines_highest_ypos(:,2));
    distSortedHlines_comb = [distSortedHlines,distSortedHlines_idx];
    diffyDistHlines = diff(distSortedHlines_comb(:,1));
    [diffyDistHlinesMax,diffyDistHlinesMax_idx] = max(diffyDistHlines);

    midptMaxDiff_Hlines_highest_ypos = (distSortedHlines(diffyDistHlinesMax_idx+1) + distSortedHlines(diffyDistHlinesMax_idx))/2;

    % Grouping
    Hlines_highest_grouped = NaN(Hlines_highest_num,1);

    for group_idx = 1:Hlines_highest_num
        if [Hlines_highest(group_idx).point1(1,2)] > midptMaxDiff_Hlines_highest_ypos
            Hlines_highest_grouped(group_idx) = 1;
            Hlines_highest(group_idx).group = 1;
        elseif [Hlines_highest(group_idx).point1(1,2)] < midptMaxDiff_Hlines_highest_ypos
            Hlines_highest_grouped(group_idx) = 2;
            Hlines_highest(group_idx).group = 2;
        else
            Hlines_highest_grouped = NaN; %error, should never happen
        end
    end

    % Group count
    [B,BG] = groupcounts(Hlines_highest_grouped);

    %% 11) Grouping Houghlines [OLD]

%     % Calculating mean ypos
%     Hlines_highest_ypos = NaN(Hlines_highest_num,3);
% 
%     for Hlines_highest_ypos_idx = 1:Hlines_highest_num
%         Hlines_highest_ypos(Hlines_highest_ypos_idx,1) = Hlines_highest_ypos_idx;
%         Hlines_highest_ypos(Hlines_highest_ypos_idx,2) = [Hlines_highest(Hlines_highest_ypos_idx).point1(1,2)];
%         Hlines_highest_ypos(Hlines_highest_ypos_idx,3) = [Hlines_highest(Hlines_highest_ypos_idx).votes];
%     end
% 
%     mean_Hlines_highest_ypos = mean(Hlines_highest_ypos(:,2));
% 
%     % Grouping
%     Hlines_highest_grouped = NaN(Hlines_highest_num,1);
% 
% %     for group_idx = 1:Hlines_highest_num
% %         if [Hlines_highest(group_idx).point1(1,2)] > mean_Hlines_highest_ypos
% %             Hlines_highest_grouped(group_idx) = 1;
% %             Hlines_highest(group_idx).group = 1;
% %         elseif [Hlines_highest(group_idx).point1(1,2)] < mean_Hlines_highest_ypos
% %             Hlines_highest_grouped(group_idx) = 2;
% %             Hlines_highest(group_idx).group = 2;
% %         else
% %             Hlines_highest_grouped = NaN; %error
% %         end
% %     end
% 
%     for group_idx = 1:Hlines_highest_num
%         if [Hlines_highest(group_idx).point1(1,2)] >= mean_Hlines_highest_ypos
%             Hlines_highest_grouped(group_idx) = 1;
%             Hlines_highest(group_idx).group = 1;
%         else
%             Hlines_highest_grouped(group_idx) = 2;
%             Hlines_highest(group_idx).group = 2;
%         end
%     end
% 
%     % Group count
%     [B,BG] = groupcounts(Hlines_highest_grouped);

    %% Houghline filtering
    
    if BG(1) == 1
        B_top_idx = 1;
        B_bot_idx = 2;
    else
        B_top_idx = 2;
        B_bot_idx = 1;
    end

    if B(B_top_idx) < 1
        warning('B2K - No top houghline detected')
        HLcorrectionidx = HLcorrectionidx + 1;
        continue
    elseif B(B_bot_idx) < 1
        warning('B2K - No bottom houghline detected')
        HLcorrectionidx = HLcorrectionidx + 1;
        continue
    elseif B(B_top_idx) < 1 || B(B_bot_idx) < 1
        warning('B2K - No houghlines detected')
        HLcorrectionidx = HLcorrectionidx + 1;
        continue
    else
        break
    end
    
end


    %% Filter for proper angles only




%%
varargout{1} = resizeScale; 
varargout{2} = edgeMethod;
varargout{3} = R_resolution; 
varargout{4} = T_resolution; 
varargout{5} = H_numpeaks;
varargout{6} = Hlines_fillgap;
varargout{7} = Hlines_minlength;
varargout{8} = Hlines_highest;
varargout{9} = Hlines_highest_num;

%% Initial figures
% Notes: tiledlayout with imref2d changes the size of the last tile if
% there is insufficient space. 

myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);

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
pause(5)
myFigures.fig.(myFigures.handle(myFigures.figNum)).WindowState ='maximized';
pause(5)

tile1 = tiledlayout(myFigures.fig.(myFigures.handle(myFigures.figNum)),14,1,'TileSpacing','tight','padding','tight');
tile1.TileIndexing = 'rowmajor';
tile1.Units = 'normalized';
tile1.OuterPosition = [0 0 1 1];
Title1 = title(tile1,'Image calibration','Interpreter','none');
Title1.FontSize = 30;

%% 1) I [TILE1]

ax1 = nexttile;
img1 = imshow(I);
img1size = size(img1.CData);
txt1 = title('I');
fontsize(txt1,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img1size(1)/2)],'HorizontalAlignment','Right');

%% 2) Edge detection, I [TILE2]

ax2 = nexttile;
img2 = imshow(edge_I);
img2size = size(img2.CData);
txt2 = title('Sobel Edge, I');
fontsize(txt2,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img2size(1)/2)],'HorizontalAlignment','Right');

%% 3) I_resize [TILE3]

ax3 = nexttile;
img3 = imshow(I_resize);
img3size = size(img3.CData);
txt3 = title('I_{resize}');
fontsize(txt3,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img3size(1)/2)],'HorizontalAlignment','Right');

%% 4) Edge detection, I_resize [TILE4]

ax4 = nexttile;
img4 = imshow(edge_I_resize);
img4size = size(img4.CData);
txt4 = title('Sobel Edge, I_{resize}');
fontsize(txt4,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img4size(1)/2)],'HorizontalAlignment','Right');

%% 5) I_average [TILE5]

ax5 = nexttile;
img5 = imshow(I_average);
img5size = size(img5.CData);
txt5 = title('I_{average}');
fontsize(txt5,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img5size(1)/2)],'HorizontalAlignment','Right');

%% 6) Edge detection, I_average [TILE6]

ax6 = nexttile;
img6 = imshow(edge_I_average);
img6size = size(img6.CData);
txt6 = title('Sobel Edge, I_{average}');
fontsize(txt6,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img6size(1)/2)],'HorizontalAlignment','Right');

%% 7) Houghlines, highest [TILE7]

ax7 = nexttile;
img7 = imshow(I_average);
img7size = size(img7.CData);
txt7 = title('houghlines,highest');
fontsize(txt7,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img7size(1)/2)],'HorizontalAlignment','Right');
hold on;
for k = 1:length(Hlines_highest)
    Hlines_full_xy_px = [Hlines_highest(k).point1; Hlines_highest(k).point2];
    plot(Hlines_full_xy_px(:,1),Hlines_full_xy_px(:,2),'LineWidth',2,'Color','green');

    % Plot beginnings and ends of lines
    plot(Hlines_full_xy_px(1,1),Hlines_full_xy_px(1,2),'x','LineWidth',2,'Color','yellow');
    plot(Hlines_full_xy_px(2,1),Hlines_full_xy_px(2,2),'x','LineWidth',2,'Color','red');
end
hold off;


%% 8) I_rotcrop [TILE8]

ax8 = nexttile;
img8 = imshow(I_rotcrop);
img8size = size(img8.CData);
txt8 = title('I_{rotcrop}');
fontsize(txt8,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img8size(1)/2)],'HorizontalAlignment','Right');

%% 9) Edge detection, I_rotcrop [TILE9]

ax9 = nexttile;
img9 = imshow(edge_I_rotcrop);
img9size = size(img9.CData);
txt9 = title('I_{rotcrop}');
fontsize(txt9,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img9size(1)/2)],'HorizontalAlignment','Right');

%% 10) Houghlines, I_rotcrop [TILE10]

ax10 = nexttile;
img10 = imshow(I_average);
img10size = size(img10.CData);
txt10 = title('houghlines,highest');
fontsize(txt10,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img10size(1)/2)],'HorizontalAlignment','Right');
hold on;
for k = 1:length(Hlines_highest)
    Hlines_full_xy_px = [Hlines_highest(k).point1; Hlines_highest(k).point2];
    plot(Hlines_full_xy_px(:,1),Hlines_full_xy_px(:,2),'LineWidth',2,'Color','green');

    % Plot beginnings and ends of lines
    plot(Hlines_full_xy_px(1,1),Hlines_full_xy_px(1,2),'x','LineWidth',2,'Color','yellow');
    plot(Hlines_full_xy_px(2,1),Hlines_full_xy_px(2,2),'x','LineWidth',2,'Color','red');
end
hold off;

%%
top_first_idx = find([Hlines_highest.group]==1,1,'first');
bot_first_idx = find([Hlines_highest.group]==2,1,'first');
t1 = 2*top_first_idx-1;
t2 = 2*top_first_idx;
b1 = 2*bot_first_idx-1;
b2 = 2*bot_first_idx;


%% 11) Intrinsic (px)

ax11 = nexttile;
img11 = imshow(I_average);
img11size = size(img11.CData);
hold on;

hough_x_px = zeros(1,2*length(Hlines_highest));
hough_y_px = zeros(1,2*length(Hlines_highest));

j = 1;
for k = 1:length(Hlines_highest)
    hough_x_px(j:j+1) = [Hlines_highest(k).point1(1),Hlines_highest(k).point2(1)];
    hough_y_px(j:j+1) = [Hlines_highest(k).point1(2),Hlines_highest(k).point2(2)];
    j=j+2;
end

plot(hough_x_px(t1:t2),hough_y_px(t1:t2),'LineWidth',2,'Color','green');
plot(hough_x_px(b1:b2),hough_y_px(b1:b2),'LineWidth',2,'Color','green');

% Plot beginnings and ends of lines
plot(hough_x_px(t1),hough_y_px(t1),'x','LineWidth',2,'Color','yellow');
plot(hough_x_px(b1),hough_y_px(b1),'x','LineWidth',2,'Color','yellow');
plot(hough_x_px(t2),hough_y_px(t2),'x','LineWidth',2,'Color','red');
plot(hough_x_px(b2),hough_y_px(b2),'x','LineWidth',2,'Color','red');

c_p1_x_px = (hough_x_px(t1) + hough_x_px(b1))/2;
c_p1_y_px = (hough_y_px(t1) + hough_y_px(b1))/2;
c_p2_x_px = (hough_x_px(t2) + hough_x_px(b2))/2;
c_p2_y_px = (hough_y_px(t2) + hough_y_px(b2))/2;

% centerline = line([c_p1_x_px c_p2_x_px],[c_p1_y_px c_p2_y_px]);
centerline = line([1-0.5 size(I_average,2)-0.5],[c_p1_y_px c_p2_y_px]);
centerline.Color = 'green';
centerline.LineStyle = '--';

axis on;
ax = gca;
% ax.XTick = 0:50:700;

ax.FontSize = fontSizeYLabel;

txt = title('intrinsic [px]');
fontsize(txt,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img11size(1)/2)],'HorizontalAlignment','Right');

hold off;

%% 12) World (mm)

p1 = [hough_x_px(t1),hough_y_px(t1),0];
p2 = [hough_x_px(t2),hough_y_px(t2),0];
p3 = [hough_x_px(b1),hough_y_px(b1),0];
p4 = [hough_x_px(b2),hough_y_px(b2),0];

[dist, vector, point1, point2] = B2KDistbw2Lines(p1, p2, p3, p4);

%known channel width for mm/px conversion
mm_per_px = EXP_W_channel/dist;

WR = imref2d(size(I_average),mm_per_px,mm_per_px);

ax12 = nexttile;
img12 = imshow(I_average,WR);
img12size = size(img12.CData);
hold on;

hough_x_px = zeros(1,length(Hlines_highest));
hough_y_px = zeros(1,length(Hlines_highest));

j = 1;
for k = 1:length(Hlines_highest)
    hough_x_px(j:j+1) = [Hlines_highest(k).point1(1),Hlines_highest(k).point2(1)];
    hough_y_px(j:j+1) = [Hlines_highest(k).point1(2),Hlines_highest(k).point2(2)];
    j=j+2;
end

[hough_x_intrinsic,hough_y_intrinsic] = intrinsicToWorld(WR,hough_x_px,hough_y_px);


plot(hough_x_intrinsic(t1:t2),hough_y_intrinsic(t1:t2),'LineWidth',2,'Color','green');
plot(hough_x_intrinsic(b1:b2),hough_y_intrinsic(b1:b2),'LineWidth',2,'Color','green');

% Plot beginnings and ends of lines
plot(hough_x_intrinsic(t1),hough_y_intrinsic(t1),'x','LineWidth',2,'Color','yellow');
plot(hough_x_intrinsic(b1),hough_y_intrinsic(b1),'x','LineWidth',2,'Color','yellow');
plot(hough_x_intrinsic(t2),hough_y_intrinsic(t2),'x','LineWidth',2,'Color','red');
plot(hough_x_intrinsic(b2),hough_y_intrinsic(b2),'x','LineWidth',2,'Color','red');

c_p1_x_intrinsic = (hough_x_intrinsic(t1) + hough_x_intrinsic(b1))/2;
c_p1_y_intrinsic = (hough_y_intrinsic(t1) + hough_y_intrinsic(b1))/2;
c_p2_x_intrinsic = (hough_x_intrinsic(t2) + hough_x_intrinsic(b2))/2;
c_p2_y_intrinsic = (hough_y_intrinsic(t2) + hough_y_intrinsic(b2))/2;


% centerline = line([c_p1_x_intrinsic c_p2_x_intrinsic],[c_p1_y_intrinsic c_p2_y_intrinsic]);
centerline = line([1-0.5 size(I_average,2)-0.5],[c_p1_y_intrinsic c_p2_y_intrinsic]);

centerline.Color = 'green';
centerline.LineStyle = '--';

axis on;
ax = gca;
% ax.XTick = 0:5:50;

% plot(c_p1_x_px,c_p1_y_px,c_p2_x_px,c_p2_y_px,'x','LineWidth',2,'Color','green');
%centerline = insertShape(img,'Line',[c_p1_x_px,c_p1_y_px,c_p2_x_px,c_p2_y_px]);

ax.FontSize = fontSizeYLabel;

hold off;

txt = title('world [mm]');
fontsize(txt,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img12size(1)/2)],'HorizontalAlignment','Right');

%% 13) Centered

CR_xWorldLimits = [WR.XWorldLimits(1)-WR.XWorldLimits(1) WR.XWorldLimits(2)-WR.XWorldLimits(1)];
CR_yWorldLimits = [WR.YWorldLimits(1)-c_p2_y_intrinsic WR.YWorldLimits(2)-c_p2_y_intrinsic];

CR = imref2d(size(I_average),CR_xWorldLimits,CR_yWorldLimits);

ax13 = nexttile;
img13 = imshow(I_average,CR);
img13size = size(img13.CData);
hold on;

hough_x_px = zeros(1,length(Hlines_highest));
hough_y_px = zeros(1,length(Hlines_highest));

j = 1;
for k = 1:length(Hlines_highest)
    hough_x_px(j:j+1) = [Hlines_highest(k).point1(1),Hlines_highest(k).point2(1)];
    hough_y_px(j:j+1) = [Hlines_highest(k).point1(2),Hlines_highest(k).point2(2)];
    j=j+2;
end

[hough_x_intrinsic,hough_y_intrinsic] = intrinsicToWorld(CR,hough_x_px,hough_y_px);


plot(hough_x_intrinsic(t1:t2),hough_y_intrinsic(t1:t2),'LineWidth',2,'Color','green');
plot(hough_x_intrinsic(b1:b2),hough_y_intrinsic(b1:b2),'LineWidth',2,'Color','green');

% Plot beginnings and ends of lines
plot(hough_x_intrinsic(t1),hough_y_intrinsic(t1),'x','LineWidth',2,'Color','yellow');
plot(hough_x_intrinsic(b1),hough_y_intrinsic(b1),'x','LineWidth',2,'Color','yellow');
plot(hough_x_intrinsic(t2),hough_y_intrinsic(t2),'x','LineWidth',2,'Color','red');
plot(hough_x_intrinsic(b2),hough_y_intrinsic(b2),'x','LineWidth',2,'Color','red');

c_p1_x_intrinsic = (hough_x_intrinsic(t1) + hough_x_intrinsic(b1))/2;
c_p1_y_intrinsic = (hough_y_intrinsic(t1) + hough_y_intrinsic(b1))/2;
c_p2_x_intrinsic = (hough_x_intrinsic(b1) + hough_x_intrinsic(b2))/2;
c_p2_y_intrinsic = (hough_y_intrinsic(b1) + hough_y_intrinsic(b2))/2;


centerline = line([1-0.5 size(I_average,2)-0.5],[c_p1_y_intrinsic c_p2_y_intrinsic]);
centerline.Color = 'green';
centerline.LineStyle = '--';

axis on;
ax = gca;
% ax.XTick = 0:5:50;

% plot(c_p1_x_px,c_p1_y_px,c_p2_x_px,c_p2_y_px,'x','LineWidth',2,'Color','green');
%centerline = insertShape(img,'Line',[c_p1_x_px,c_p1_y_px,c_p2_x_px,c_p2_y_px]);

ax.FontSize = fontSizeYLabel;

hold off;
txt = title('centered [mm]');
fontsize(txt,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img13size(1)/2)],'HorizontalAlignment','Right');

%% 14) Normalized

%known width channel for conversion
norm_per_mm = 1/(EXP_W_channel/2);

NR_xWorldLimits = [CR.XWorldLimits(1) CR.XWorldLimits(2)];
NR_yWorldLimits = [CR.YWorldLimits(1)*(norm_per_mm) CR.YWorldLimits(2)*(norm_per_mm)];

% Tried to make height uniform, but the scaling yieled non-integer values
% sizeI = size(I);
% sizeNRy = sizeI(1)*norm_per_mm;
% sizeNRx = sizeI(2);
% sizeNR = [sizeNRy sizeNRx];
% 
% NR = imref2d(sizeNR,NR_xWorldLimits,NR_yWorldLimits);

NR = imref2d(size(I_average),NR_xWorldLimits,NR_yWorldLimits);

ax14 = nexttile;
img14 = imshow(I_average,NR);
img14size = size(img13.CData);
hold on;

hough_x_px = zeros(1,length(Hlines_highest));
hough_y_px = zeros(1,length(Hlines_highest));

j = 1;
for k = 1:length(Hlines_highest)
    hough_x_px(j:j+1) = [Hlines_highest(k).point1(1),Hlines_highest(k).point2(1)];
    hough_y_px(j:j+1) = [Hlines_highest(k).point1(2),Hlines_highest(k).point2(2)];
    j=j+2;
end

[hough_x_intrinsic,hough_y_intrinsic] = intrinsicToWorld(NR,hough_x_px,hough_y_px);


plot(hough_x_intrinsic(t1:t2),hough_y_intrinsic(t1:t2),'LineWidth',2,'Color','green');
plot(hough_x_intrinsic(b1:b2),hough_y_intrinsic(b1:b2),'LineWidth',2,'Color','green');

% Plot beginnings and ends of lines
plot(hough_x_intrinsic(t1),hough_y_intrinsic(t1),'x','LineWidth',2,'Color','yellow');
plot(hough_x_intrinsic(b1),hough_y_intrinsic(b1),'x','LineWidth',2,'Color','yellow');
plot(hough_x_intrinsic(t2),hough_y_intrinsic(t2),'x','LineWidth',2,'Color','red');
plot(hough_x_intrinsic(b2),hough_y_intrinsic(b2),'x','LineWidth',2,'Color','red');

c_p1_x_intrinsic = (hough_x_intrinsic(t1) + hough_x_intrinsic(b1))/2;
c_p1_y_intrinsic = (hough_y_intrinsic(t1) + hough_y_intrinsic(b1))/2;
c_p2_x_intrinsic = (hough_x_intrinsic(b1) + hough_x_intrinsic(b2))/2;
c_p2_y_intrinsic = (hough_y_intrinsic(b1) + hough_y_intrinsic(b2))/2;


centerline = line([1-0.5 size(I_average,2)-0.5],[c_p1_y_intrinsic c_p2_y_intrinsic]);
centerline.Color = 'green';
centerline.LineStyle = '--';

axis on;
ax = gca;
% ax.XTick = 0:5:50;
ax.YTick = -1:1:1;

% plot(c_p1_x_px,c_p1_y_px,c_p2_x_px,c_p2_y_px,'x','LineWidth',2,'Color','green');
%centerline = insertShape(img,'Line',[c_p1_x_px,c_p1_y_px,c_p2_x_px,c_p2_y_px]);

ax.FontSize = fontSizeYLabel;

hold off;
txt = title('normalized [-]');
fontsize(txt,txtFontSizePixel,"pixels");
set(get(gca,'title'),'Units','Normalized','Position',[tile_x_pos_norm (0.5-txtFontSizePixel/img14size(1)/2)],'HorizontalAlignment','Right');

%% Coordinate System
% IR, WR, CR, NR

switch CS_string
    case 'IR' %Intrinsic
        x_pos_caption = ' X-pos (px):';
        y_pos_caption = ' Y-pos (px):';
        CS = IR;
    case 'WR' %World
        x_pos_caption = ' X-pos (mm):';
        y_pos_caption = ' Y-pos (mm):';
        CS = WR;
    case 'CR' %Centered
        x_pos_caption = ' X-pos (mm):';
        y_pos_caption = ' Y-pos (mm):';
        CS = CR;
    case 'NR' %Normalized
        x_pos_caption = ' X-pos (mm):';
        y_pos_caption = ' Y-pos (-):';
        CS = NR;
end

%%
clear B2KMyFigures

end
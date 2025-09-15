function [myFigures,position_x_CS,position_y_CS] = B2KParticleTracking03_orig(i,vidObj,myFigures,onePTrackOutputFolderPath,L_start,R_end,tFrames,nFrames,I,CS,x_pos_caption,y_pos_caption,EXP_actualframerate,flagDebug,flagRunWithDisplay)
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

% %% Batch data table setup
% batchdata = table();
% 
% %% Data export setup
% outputfolder = 'E:\OutputFolder';
% 
% %% Video setup
% 
% % vid = '1p_Isopropanol_R1_221109_H1.0_W1.5_4mL-min_1500FPS_1.avi';
% 
% vid = '1p_Water_R1_221026_H1.0_W1.5_25mL-min_1500FPS_1.avi';
% 
% % vid = 'Water_R1_230125_H1.0_W1.5_15mL-min_1500FPS_9.avi';
% 
% % vid = 'Water_R1_221101_H1.5_W1.0_20mL-min_1500FPS_1.avi'; %[problem]
% % vid = 'Water_R1_221026_H1.0_W1.5_9mL-min_1500FPS_2.avi';
% vidObj = VideoReader(vid);

%% Video path, video writer object setup
vidNumStr = sprintf('_vid%d',i);
vidOutputOrigFileName = 'B2K_p-Chip_Tracking.avi';
vidOutputFileName = insertBefore(vidOutputOrigFileName,'.avi',vidNumStr);
file = fullfile(onePTrackOutputFolderPath,vidOutputFileName);
vidWrtObj = VideoWriter(file,'Uncompressed AVI');
vidWrtObj.FrameRate = 15;
open(vidWrtObj);

%% Tracking setup
L_start_orig = L_start;
R_end_orig = R_end;

%% Plot variables
tile_x_pos = -200;
tile_y_pos = 30;

%% Create vision.BlobAnalysis system object
% sysObj_Vision_BlobAnalysis = vision.BlobAnalysis('MinimumBlobArea',100,'MaximumBlobArea',500,'PerimeterOutputPort',true,'OrientationOutputPort',true,'ExcludeBorderBlobs',true);
sysObj_Vision_BlobAnalysis = vision.BlobAnalysis('MinimumBlobArea',50,'MaximumBlobArea',500,'PerimeterOutputPort',true,'OrientationOutputPort',true,'ExcludeBorderBlobs',true); %error with B2KMinBoundParalellogram
% sysObjVision_BlobAnalysis = vision.BlobAnalysis('MinimumBlobArea',75,'MaximumBlobArea',500,'PerimeterOutputPort',true,'OrientationOutputPort',true,'ExcludeBorderBlobs',true); %error with B2KMinBoundParalellogram

%% Preallocate all frames for image binarization, connected components, regionprops
binaryframes = zeros(size(I,1),size(I,2),nFrames);
CC_store = cell(nFrames,1);
S_store = cell(nFrames,1);

%% Preallocate for position, velocity, acceleration data for each frames
position_x_CS = NaN([1,tFrames],'double');
velocity_x_CS = NaN([1,tFrames],'double');
acceleration_x_CS = NaN([1,tFrames],'double');

position_y_CS = NaN([1,tFrames],'double');
velocity_y_CS = NaN([1,tFrames],'double');
acceleration_y_CS = NaN([1,tFrames],'double');

velocity_y_CS_avg = NaN([1,tFrames],'double');
acceleration_y_CS_avg = NaN([1,tFrames],'double');

frame_array = 1:1:tFrames;
frame_num_array = L_start_orig:R_end_orig;
time_ms_array = frame_array*(1/EXP_actualframerate)*(1000);

%% Preallocate blobCentroid and Harris Centroid data (comparing both methods for centroid detection)
blobCentroid_x = NaN([1,tFrames],'double');
blobCentroid_y = NaN([1,tFrames],'double');
HarrisCentroid_x = NaN([1,tFrames],'double');
HarrisCentroid_y = NaN([1,tFrames],'double');

blobCentroid_x_CS = NaN([1,tFrames],'double');
blobCentroid_y_CS = NaN([1,tFrames],'double');
HarrisCentroid_x_CS = NaN([1,tFrames],'double');
HarrisCentroid_y_CS = NaN([1,tFrames],'double');

centroid_distance_CS = NaN([1,tFrames],'double');

%% Preallocate objTracks/cornerTracks and cell for storing objTracks/cornerTracks
objTracks = table(); %Creating table for objTracks; enables loop to run when there are no objTracks yet
cornerTracks = table(); %Creating table for cornerTracks; enables loop to run when there are no cornerTracks yet

storedObjTracks = cell(tFrames,1);
storedCornerTracks = cell(tFrames,1);

%
myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);

%L_start_orig = 1269;
%L_firstDetection = 1307; %row 39 in storedObjTracks
%L_firstLostDetection = 1345; (found again 1346)
%L_secondLostDetection = 1368; (found again 1369, lost 1370, found again 1399)
% Loop through frames between the start frame and end frame evaluated from Video Summarization
for i = L_start:R_end 
    frame = read(vidObj,i); %Read designated frame from selected video
    gray = rgb2gray(frame); %Convert RGB video (3D matrix) to grayscale (2D matrix); required even though video is in grayscale since the video is saved as 3D matrix
    W_gray = wiener2(gray,[3,3]); %Filter using Wiener filtering
    Lap_W_gray = locallapfilt(W_gray,0.09,5); %Additional filter using Laplacian filtering
    binaryframes(:,:,i) = ~imbinarize(Lap_W_gray,0.32); %Image binarization
    CC = bwconncomp(binaryframes(:,:,i), 8); %Find & count connected components in the binary image
    CC_store(i,1) = {CC};
    S = regionprops(CC,'all');
    S_store(i,1) = {S};
    binaryframes(:,:,i) = bwareaopen(binaryframes(:,:,i),50); %Remove small objects from binary image
    % binaryframes_L = logical(~binaryframes);
    % binaryframes_L = ~binaryframes;
    % binaryframes_L = imcomplement(binaryframes);
    binaryframes_L = ~binaryframes; %Change object(s) to be the white pixels and the background to black pixels

    [objArea,objCentroid,bboxOut,perimeter,orientation] = step(sysObj_Vision_BlobAnalysis, ~binaryframes_L(:,:,i)); %Perform Blob Analysis to find "blobs" corresponding to particle(s) (& inadvertently artifacts)

    numObj = length(objArea); %blob number

    bboxOut_spacing = [-5 -5 10 10];

    % If there is a detected blob object, then:
    if isempty(bboxOut)
        largerbboxOut = zeros(1,4); %Create a dummy object for a non-existent bounding box
        Ishape2 = insertShape(Lap_W_gray,'rectangle',largerbboxOut,'Linewidth',1,'Color','yellow','Opacity',0.001); %Set the finalized shape with a non-existent bounding box; this is required 
        % Ishape2 = Lap_W_gray; %this replacement line runs properly and saves all output images/videos, but hits a closing figure handle error: Error using close Invalid figure handle. Error in Batch1pTracking (line 987): close(myFigures.fig.(myFigures.handle(idxFigRef)))

        %%
        % objDetections [1x3 table]: objCentroid [1x2: x,y], bboxOut [1x4: x,y,w,h based on top left corner], assigned [1x1: 0 or 1]
        % objTracks [1x11 table]: objTrackID [1x1: #], DetectedLocation
        % [1x2: x,y], PredictedLocation [1x2: x,y], TrackedLocation [1x2:
        % x,y], PredictedState [1x4: ????], KalmanFilter [1x1: KalmanFilter
        % object], TotalFrameCount [1x1: # of frames that object is
        % detected], TotalVisibleCount [1x1: #?], ConsecutiveInvisibleCount
        % [1x1: #?], Visible [1x1: #], Confirmed [1x1: #]

        % Notes: - need a separate conditional statement to handle initial
        % frames when the object hasn't been detected at all yet
        % - need a flag/toggle
        
    else
        largerbboxOut = bboxOut + int32(bboxOut_spacing);

        % Find artifacts in first frame
        % Binarize the image by thresholding.
        mask = ~imbinarize(Lap_W_gray,0.38);
        % % Find the areas
        % props = regionprops(mask, 'Area');
        % allAreas = sort([props.Area])

        % % Extract only blobs larger than 100.
        % mask = bwareafilt(mask,[5 100],8);
        % Lap_W_gray(mask) = mean(Lap_W_gray(1,:)); % Now mask

        % Extract only blobs larger than 50.
        mask = bwareafilt(mask,[5 50],8);
        Lap_W_gray(mask) = mean(Lap_W_gray(1,:)); % Now mask

        % subpixel calculation
        threshold = 10;
        edges = subpixelEdges(Lap_W_gray, threshold);
        ex = edges.x;
        ey = edges.y;
        metric = 'a';
        [pgx,pgy,area,perim] = B2KMinBoundParallelogram(ex,ey,metric);
        % [pgx,pgy,~,~] = B2KMinBoundParallelogram(ex,ey,metric);
        Ishape = insertShape(Lap_W_gray,'polygon',[pgx,pgy],'Linewidth',1,'Color','y','Opacity',0.001);
        Ishape2 = insertShape(Ishape,'rectangle',largerbboxOut,'Linewidth',1,'Color','g','Opacity',0.001);
        
        %%
        % Object Detection 
        objDetections = detectObjs(objCentroid,bboxOut); %Immediately outputs values for objDetections

        % % Store previous object detection
        % prevobjDetections = storeprevobjDetections(objDetections);

        % Object Prediction [last known position]
        objTracks = predictobjTracks(objTracks); %Immediately updates PredictedLocation in objTracks

        % Object Assign Tracks to Detections
        [objTracks, objDetections] = assignTracksDetections(objDetections,objTracks); %Immediately updates DetectedLocation in objTracks

        % Update objTracks
        objTracks = updateobjTracks(objTracks,objDetections); %Immediately updates TrackedLocation in objTracks
        
        %%
        % Corner Detection
        cornerDetections = detectCorners(pgx,pgy);

        % Corner Prediction
        cornerTracks = predictcornerTracks(cornerTracks);

        % Corner Assign Object Tracks & Corner Tracks to Corner Detections
        [cornerTracks, cornerDetections] = assigncornerTracksDetections(cornerTracks,cornerDetections);

        % Update cornerTracks
        cornerTracks = updatecornerTracks(cornerTracks,cornerDetections);

    end

    text0origin = [0 60];
    text1origin = [0 0];
    text2spacing = [80 0];
    text2addspacing = [80 0];

    Itext0 = insertText(Ishape2,text0origin,['Captured FPS: ' num2str(EXP_actualframerate) '   |   Playback speed: ' num2str(EXP_actualframerate/vidWrtObj.FrameRate,'%.2f') 'X slower'],'FontSize',10,'BoxColor','red','BoxOpacity',0.1,'TextColor','white'); %bottom left red text box describing capture FPS & playback in exported video
    Itext1 = insertText(Itext0,text1origin,['Frame #: ' num2str(i)],'FontSize',10,'BoxColor','green','BoxOpacity',0.1,'TextColor','white'); %top left green text box describing Frame # in exported video

    if isempty(objCentroid)
        objCentroid = [0 0];
        Itext2 = insertText(Itext1,[80 0],'No p-Chip found','FontSize',10,'BoxColor','blue','BoxOpacity',0.2,'TextColor','white'); %top left blue text box describing p-Chip detection/tracking information in exported video
    else %if numObj >=1

        % Transform data to proper coordinate system
        [objCentroid_x_CS,objCentroid_y_CS] = intrinsicToWorld(CS,objCentroid(1,1),objCentroid(1,2));

        % Store objCentroid data in blobCentroid array
        blobCentroid_x_CS(i - L_start_orig + 1) = objCentroid_x_CS;
        blobCentroid_y_CS(i - L_start_orig + 1) = objCentroid_y_CS;

        % Position calculation [mm]
        position_x_CS(i - L_start_orig + 1) = objCentroid_x_CS;
        position_y_CS(i - L_start_orig + 1) = objCentroid_y_CS;

        firstframe = frame_num_array(find(~isnan(position_x_CS),1));

        if i == firstframe
            velocity_x_CS(i - L_start_orig + 1) = NaN;
            acceleration_x_CS(i - L_start_orig + 1) = NaN;
            velocity_y_CS(i - L_start_orig + 1) = NaN;
            acceleration_y_CS(i - L_start_orig + 1) = NaN;
        elseif i == firstframe + 1
            velocity_x_CS(i - L_start_orig + 1) = (position_x_CS(i - L_start_orig + 1) - position_x_CS(i - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]
            acceleration_x_CS(i - L_start_orig + 1) = NaN;
            velocity_y_CS(i - L_start_orig + 1) = (position_y_CS(i - L_start_orig + 1) - position_y_CS(i - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]
            acceleration_y_CS(i - L_start_orig + 1) = NaN;
        else
            % Velocity calculation [m/s]
            velocity_x_CS(i - L_start_orig + 1) = (position_x_CS(i - L_start_orig + 1) - position_x_CS(i - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]
            velocity_y_CS(i - L_start_orig + 1) = (position_y_CS(i - L_start_orig + 1) - position_y_CS(i - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]

            % Acceleration calculation [m/s^2]
            acceleration_x_CS(i - L_start_orig + 1) = (velocity_x_CS(i - L_start_orig + 1) - velocity_x_CS(i - L_start_orig))*(EXP_actualframerate); %[m/s^2]
            acceleration_y_CS(i - L_start_orig + 1) = (velocity_y_CS(i - L_start_orig + 1) - velocity_y_CS(i - L_start_orig))*(EXP_actualframerate); %[m/s^2]
        end

        %Note: Take movmean for velocity either after the end or do the movmean at real time 
        % Multi-point average of velocity on x & y by 5 points
        velocity_x_CS_avg = movmean(velocity_x_CS, 5);
        velocity_y_CS_avg = movmean(velocity_y_CS, 5);

        % Multi-point average of acceleration on x & y by 5 points
        acceleration_x_CS_avg = movmean(acceleration_x_CS, 5);
        acceleration_y_CS_avg = movmean(acceleration_y_CS, 5);

        % Calculate and store HarrisCentroid data in array
        points = detectHarrisFeatures(binaryframes_L(:,:,i));
        corners = points.selectStrongest(4);
        % Transform data to proper coordinate system
        [corners_x_CS,corners_y_CS] = intrinsicToWorld(CS,corners.Location(:,1),corners.Location(:,2));

        % Insert color-coded corner markers
        cornermarker1 = insertMarker(Itext1,cornerTracks.DetectedLocation(1,:),'+','Color','red');
        cornermarker2 = insertMarker(cornermarker1,cornerTracks.DetectedLocation(2,:),'+','Color','green');
        cornermarker3 = insertMarker(cornermarker2,cornerTracks.DetectedLocation(3,:),'+','Color','blue');
        cornersmarkers = insertMarker(cornermarker3,cornerTracks.DetectedLocation(4,:),'+','Color','magenta');

        % Determine Harris corner centroid
        HarrisCentroid_x_CS(i - L_start_orig + 1) = mean(corners_x_CS);
        HarrisCentroid_y_CS(i - L_start_orig + 1) = mean(corners_y_CS);
        HarrisCentroid_x(i - L_start_orig + 1) = mean(corners.Location(:,1));
        HarrisCentroid_y(i - L_start_orig + 1) = mean(corners.Location(:,2));

        % Calculate distance between blobCentroid and HarrisCentroid
        centroid_distance_CS = sqrt((blobCentroid_x_CS - HarrisCentroid_x_CS).^2+(blobCentroid_y_CS - HarrisCentroid_y_CS).^2);

        text_str = cell(numObj,1);
        text2origin = zeros(numObj,2);

        % [Old iterative figure plotting]
        for ii = 1:numObj
            text_str{ii} = ['[p-Chip' num2str(ii) ']' x_pos_caption num2str(objCentroid_x_CS,'%.2f') ',' y_pos_caption num2str(objCentroid_y_CS,'%.2f')];
            text2origin(ii,:) = text1origin + text2spacing;
            text2spacing = text2spacing + text2addspacing;
        end
        Itext2 = insertText(cornersmarkers,text2origin,text_str,'FontSize',10,'BoxColor','blue','BoxOpacity',0.2,'TextColor','white');
    end

    %Display image
    if flagRunWithDisplay == 1
        imshow(Itext2,'InitialMagnification',500); %Can suppress this to improve performance
    end
    writeVideo(vidWrtObj,Itext2);

    % Store
    storedObjTracks{i-L_start_orig+1,1} = objTracks;
    storedCornerTracks{i-L_start_orig+1,1} = cornerTracks;

end

%%
myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);
figure(myFigures.fig.(myFigures.handle(myFigures.figNum)))
figXPosData = gcf;
tileXPosData = tiledlayout(figXPosData,3,1,'TileSpacing','tight','padding','tight');
tileXPosData.TileIndexing = 'rowmajor';
tileXPosData.OuterPosition = [0.25 0 0.5 1];
TitleXPosData = title(tileXPosData,'Position Data','Interpreter','none');
TitleXPosData.FontSize = 30;

nexttile(1)
plot(time_ms_array,position_x_CS,"o");
txt = title('X-position');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Position [mm]')
hold on

nexttile(2)
plot(time_ms_array,velocity_x_CS,"o"); %[m/s]
txt = title('X-velocity');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Velocity [m/s]')
hold on

nexttile(2)
plot(time_ms_array, velocity_x_CS_avg, "r", LineWidth=2.0); %[m/s^2]
hold on

nexttile(3)
plot(time_ms_array,acceleration_x_CS,"o"); %[m/s^2]
txt = title('X-acceleration');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Acceleration [m/s^{2}]')
hold on

nexttile(3)
plot(time_ms_array, acceleration_x_CS_avg, "r", LineWidth=2.0); %[m/s^2]
hold on

%%
myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);
figure(myFigures.fig.(myFigures.handle(myFigures.figNum)))
figYPosData = gcf;
tileYPosData = tiledlayout(figYPosData,3,1,'TileSpacing','tight','padding','tight');
tileYPosData.TileIndexing = 'rowmajor';
tileYPosData.OuterPosition = [0.25 0 0.5 1];
TitleYPosData = title(tileYPosData,'Position Data','Interpreter','none');
TitleYPosData.FontSize = 30;

nexttile(1)
plot(time_ms_array,position_y_CS,"o");
txt = title('Y-position');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Position [mm]')
hold on

nexttile(2)
plot(time_ms_array,velocity_y_CS,"o"); %[m/s]
txt = title('Y-velocity');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Velocity [m/s]')
hold on

nexttile(2)
plot(time_ms_array, velocity_y_CS_avg, "r", LineWidth=2.0); %[m/s^2]
hold on

nexttile(3)
plot(time_ms_array,acceleration_y_CS,"o"); %[m/s^2]
txt = title('Y-acceleration');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Acceleration [m/s^{2}]')
hold on

nexttile(3)
plot(time_ms_array, acceleration_y_CS_avg, "r", LineWidth=2.0); %[m/s^2]
hold on

%%
myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);
figure(myFigures.fig.(myFigures.handle(myFigures.figNum)))
figCentroidData = gcf;
tileCentroidData = tiledlayout(figCentroidData,3,1,'TileSpacing','tight','padding','tight');
tileCentroidData.TileIndexing = 'rowmajor';
tileCentroidData.OuterPosition = [0.25 0 0.5 1];
TitleCentroidData = title(tileCentroidData,'Comparison between blobCentroid and HarrisCentroid','Interpreter','none');
TitleCentroidData.FontSize = 30;

nexttile(1)
plot(L_start_orig:1:R_end_orig,HarrisCentroid_x_CS,"o");
txt = title('Centroid x-position');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('Frame # [-]')
ylabel('X-Position [mm]')
hold on
plot(L_start_orig:1:R_end_orig,blobCentroid_x_CS,"x")

nexttile(2)
plot(L_start_orig:1:R_end_orig,centroid_distance_CS,"o");
txt = title('Centroid difference');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('Frame # [-]')
ylabel('Centroid distance difference [mm]')
hold on

nexttile(3)
centroid_distance_CS_notnan = centroid_distance_CS(~isnan(centroid_distance_CS));
hist = histogram(centroid_distance_CS_notnan,"BinMethod","fd");
txt = title('Histogram of centroid dist. diff.');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('Centroid distance difference [mm]')
ylabel('Count [-]')
hold on

nexttile(2)
[counts, binEdges] = histcounts(centroid_distance_CS_notnan,"BinMethod","fd");
histIndex = find(counts < 3, 1);
centroid_threshold = binEdges(histIndex);
centroid_threshold_str = sprintf("Threshold: %.2f",centroid_threshold);
yline(centroid_threshold,'-r',centroid_threshold_str);
hold on

problemframes_idx = find(centroid_distance_CS>centroid_threshold);
problemframes = frame_num_array(problemframes_idx(1,:));

nexttile(3)
xline(centroid_threshold,'-r',centroid_threshold_str,LabelVerticalAlignment='middle');

mask = centroid_distance_CS_notnan > centroid_threshold;
histogram(centroid_distance_CS_notnan(~mask),binEdges,'FaceColor','b');
histogram(centroid_distance_CS_notnan(mask),binEdges,'FaceColor','r');

hold off

%%
myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);
figure(myFigures.fig.(myFigures.handle(myFigures.figNum)))
figDetectvsKalman = gcf;
tileDetectvsKalman = tiledlayout(figDetectvsKalman,3,1,'TileSpacing','tight','padding','tight');
tileDetectvsKalman.TileIndexing = 'rowmajor';
tileDetectvsKalman.OuterPosition = [0.25 0 0.5 1];
TitleDetectvsKalman = title(tileDetectvsKalman,'Comparison of Detections and Kalman Filter Estimates','Interpreter','none');
TitleDetectvsKalman.FontSize = 30;

nexttile
hold on;
ylabel('X-position')
xlabel('Frame')
for idx = 1:length(storedObjTracks)
    if ~isempty(storedObjTracks{idx,1})
        plot(storedObjTracks{idx,1}.TotalFrameCount,storedObjTracks{idx,1}.DetectedLocation(1,1),'k+')
        plot(storedObjTracks{idx,1}.TotalFrameCount,storedObjTracks{idx,1}.PredictedLocation(1,1),'ro')
    end
end
legend('Detected locations','Predicted/corrected locations');
hold off;

nexttile
hold on;
ylabel('Y-position')
xlabel('Frame')
for idx = 1:length(storedObjTracks)
    if ~isempty(storedObjTracks{idx,1})
        plot(storedObjTracks{idx,1}.TotalFrameCount,storedObjTracks{idx,1}.DetectedLocation(1,2),'k+')
        plot(storedObjTracks{idx,1}.TotalFrameCount,storedObjTracks{idx,1}.PredictedLocation(1,2),'ro')
    end
end
legend('Detected locations','Predicted/corrected locations');
hold off;

%% Close video writer object
close(vidWrtObj);

%% Clean up system objects
release(sysObj_Vision_BlobAnalysis)

%%
clear B2KMyFigures

%% FUNCTIONS - objDetections = detectObjs(objCentroid,bboxOut)
function objDetections = detectObjs(objCentroid,bboxOut)
    objDetections = table(objCentroid,bboxOut);
end

%% FUNCTIONS - objTracks = predictobjTracks(objTracks) [last known position]
function objTracks = predictobjTracks(objTracks)
    for idx = 1:height(objTracks)
        objTracks.PredictedLocation(idx,:) = predict(objTracks.KalmanFilter{idx}); %[1 - TRACK PREDICTION]
    end
end

%% FUNCTIONS - [objTracks, objDetections] = assignTracksDetections(objDetections,objTracks, bboxOut)
function [objTracks, objDetections] = assignTracksDetections(objDetections,objTracks, bboxOut)
    
    % Compute the cost of assigning each detection to each track [1 - ASSIGNMENT COSTS]
    cost = zeros(height(objTracks),height(objDetections));
    for idx = 1:height(objTracks)
        cost(idx,:) = distance(objTracks.KalmanFilter{idx},objDetections.objCentroid);
    end
    
    % Solve the assignment problem 
    costOfNonAssignment = 1000; %[2 - COST OF NON-ASSIGNMENT]
    assignedIdxPairs = assignDetectionsToTracks(cost,costOfNonAssignment); %[3 - ASSIGNMENT OPTIMIZATION]

    % Set objdetected flag & add objdetected location to objdetected tracks
    if ~isempty(objTracks)
        objTracks.Visible(:) = false;
        objTracks.Visible(assignedIdxPairs(:,1)) = true;
        objTracks.DetectedLocation(assignedIdxPairs(:,1),:) = objDetections.objCentroid(assignedIdxPairs(:,2),:);
    end

    % Set assignment flag to detections
    if ~isempty(objDetections)
        objDetections.assigned(:) = false;
        objDetections.assigned(assignedIdxPairs(:,2)) = true;
    end
end

%% FUNCTIONS - objTracks = updateobjTracks(objTracks,objDetections)
function objTracks = updateobjTracks(objTracks,objDetections)
    
    % Initialize persistent counter
    persistent objTrackID;
    if isempty(objTrackID)
        objTrackID = 1; %If no object, assign trackID of 1
    end

    % Update tracked locations using available detections [1 - TRACK ESTIMATES UPDATE]
    for idx = 1:height(objTracks)
        if objTracks.Visible(idx)
            % objTracks.TrackedLocation(idx,:) = predict(objTracks.KalmanFilter{idx});
            objTracks.TrackedLocation(idx,:) = correct(objTracks.KalmanFilter{idx},objTracks.DetectedLocation(idx,:)); %Provides a corrected location based on (average of?) predicted location & detected location
        end
    end

    % Update remaining track variables [1 - TRACK ESTIMATES UPDATE]
    if ~isempty(objTracks)
        % Set objTracked location to prediction for undetected objTracks
        objTracks.TrackedLocation(~objTracks.Visible,:) = objTracks.PredictedLocation(~objTracks.Visible,:);

        % Increment all track ages [2 - TRACK METADATA UPDATE]
        objTracks.TotalFrameCount = objTracks.TotalFrameCount + 1;
        
        % Increment all visible tracks total visible count [2 - TRACK METADATA UPDATE]
        objTracks.TotalVisibleCount(objTracks.Visible) = objTracks.TotalVisibleCount(objTracks.Visible) + 1;
        
        % Confirm tracks detected enough [2 - TRACK METADATA UPDATE]
        objTrackConfirmationThreshold = 3;
        objTracks.Confirmed = objTracks.TotalVisibleCount > objTrackConfirmationThreshold;

        % Reset the counter for consecutive frames undetected if detected [2 - TRACK METADATA UPDATE]
        objTracks.ConsecutiveInvisibleCount(objTracks.Visible) = 0;

        % Increment the counter for consecutive frames undetected [2 - TRACK METADATA UPDATE]
        objTracks.ConsecutiveInvisibleCount(~objTracks.Visible) = objTracks.ConsecutiveInvisibleCount(~objTracks.Visible) + 1;

        % Compute the fraction of each objTracks's age for which it was visible
        visibility = objTracks.TotalVisibleCount ./ objTracks.TotalFrameCount;
        
        % Set thresholds for age and visibility and considered missing
        ageThreshold = 10; %RY need to adjust
        visibilityThreshold = 0.6; %RY need to adjust
        lostThreshold = 10; %RY need to adjust

        % Find the indices of objTracks to delete [5 - LOST TRACK DELETION]
        newInds = objTracks.TotalFrameCount <= ageThreshold;
        lowVisibilityInds = visibility < visibilityThreshold;
        lostInds = objTracks.ConsecutiveInvisibleCount >= lostThreshold;
        
        deleteInds = (newInds & lowVisibilityInds) | lostInds;

        % Delete tracks [5 - LOST TRACK DELETION]
        objTracks = objTracks(~deleteInds,:);
    end

    % Create new objTracks for unassigned objDetections [3 - NEW TRACK INITIALIZATION]
    for idx = 1:height(objDetections)
        if ~objDetections.assigned(idx)

            % Set Kalman Filter parameters [4 - TRACK INITIALIZATION PARAMETERS]
            InitialLocation = objDetections.objCentroid(idx,:);
            FilterType = "ConstantVelocity";
            InitialEstimateError = [200,50];
            MotionNoise = [100, 25]; 
            MeasurementNoise = 100;

            % Create a Kalman Filter object [3 - NEW TRACK INITIALIZATION]
            KalmanFilter = configureKalmanFilter(FilterType,InitialLocation, ...
                InitialEstimateError, MotionNoise, MeasurementNoise);
            
            % Set initial variables for new objTrack [4 - TRACK INITIALIZATION PARAMETERS]
            DetectedLocation = InitialLocation;
            TrackedLocation = InitialLocation;
            PredictedLocation = InitialLocation;
            PredictedState = zeros(1,4);
            TotalFrameCount = 1;
            TotalVisibleCount = 1;
            ConsecutiveInvisibleCount = 0;
            Visible = true;
            Confirmed = false;

            % Convert Kalman Filter to a cell for stacking
            KalmanFilter = {KalmanFilter};

            % Create new objTrack row [4 - TRACK INITIALIZATION PARAMETERS]
            newobjTrack = table(objTrackID,...
                DetectedLocation,...
                PredictedLocation,...
                TrackedLocation,...
                PredictedState,...
                KalmanFilter,...
                TotalFrameCount,...
                TotalVisibleCount,...
                ConsecutiveInvisibleCount,...
                Visible,...
                Confirmed);
            
            % Add it to the array of tracks [4 - TRACK INITIALIZATION PARAMETERS]
            objTracks = [objTracks; newobjTrack];

            % Iterate objTrackID [4 - TRACK INITIALIZATION PARAMETERS]
            objTrackID = objTrackID + 1;
        end
    end
end

%% FUNCTIONS - objDetections = lostDetectObjs(objCentroid,bboxOut)
function objDetections = lostDetectObjs(objCentroid,bboxOut)
    objDetections = table(objCentroid,bboxOut);
end

%% FUNCTIONS - objTracks = lostPredictobjTracks(objTracks) [last known position]
function objTracks = lostPredictobjTracks(objTracks)
    for idx = 1:height(objTracks)
        objTracks.PredictedLocation(idx,:) = predict(objTracks.KalmanFilter{idx});
    end
end

%% FUNCTIONS - [objTracks, objDetections] = lostAssignTracksDetections(objDetections,objTracks)
function [objTracks, objDetections] = lostAssignTracksDetections(objDetections,objTracks)
    
    % Compute the cost of assigning each detection to each track
    cost = zeros(height(objTracks),height(objDetections));
    for idx = 1:height(objTracks)
        cost(idx,:) = distance(objTracks.KalmanFilter{idx},objDetections.objCentroid);
    end
    
    % Solve the assignment problem
    costOfNonAssignment = 1000;
    assignedIdxPairs = assignDetectionsToTracks(cost,costOfNonAssignment);

    % Set objdetected flag & add objdetected location to objdetected tracks
    if ~isempty(objTracks)
        objTracks.Visible(:) = false;
        objTracks.Visible(assignedIdxPairs(:,1)) = true;
        objTracks.DetectedLocation(assignedIdxPairs(:,1),:) = objDetections.objCentroid(assignedIdxPairs(:,2),:);
    end

    % Set assignment flag to detections
    if ~isempty(objDetections)
        objDetections.assigned(:) = false;
        objDetections.assigned(assignedIdxPairs(:,2)) = true;
    end
end

%% FUNCTIONS - objTracks = lostUpdateobjTracks(objTracks,objDetections)
function objTracks = lostUpdateobjTracks(objTracks,objDetections)
    
    % Initialize persistent counter
    persistent objTrackID;
    if isempty(objTrackID)
        objTrackID = 1;
    end

    % Update tracked locations using available detections
    for idx = 1:height(objTracks)
        if objTracks.Visible(idx)
            objTracks.TrackedLocation(idx,:) = correct(objTracks.KalmanFilter{idx},objTracks.DetectedLocation(idx,:));
        end
    end

    % Update remaining track variables
    if ~isempty(objTracks)
        % Set objTracked location to prediction for undetected objTracks
        objTracks.TrackedLocation(~objTracks.Visible,:) = objTracks.PredictedLocation(~objTracks.Visible,:);

        % Increment all track ages
        objTracks.TotalFrameCount = objTracks.TotalFrameCount + 1;
        
        % Increment all visible tracks total visible count
        objTracks.TotalVisibleCount(objTracks.Visible) = objTracks.TotalVisibleCount(objTracks.Visible) + 1;
        
        % Confirm tracks detected enough
        objTrackConfirmationThreshold = 3;
        objTracks.Confirmed = objTracks.TotalVisibleCount > objTrackConfirmationThreshold;

        % Reset the counter for consecutive frames undetected if detected
        objTracks.ConsecutiveInvisibleCount(objTracks.Visible) = 0;

        % Increment the counter for consecutive frames undetected
        objTracks.ConsecutiveInvisibleCount(~objTracks.Visible) = objTracks.ConsecutiveInvisibleCount(~objTracks.Visible) + 1;

        % Compute the fraction of each objTracks's age for which it was vis
        visibility = objTracks.TotalVisibleCount ./ objTracks.TotalFrameCount;
        
        % Set thresholds for age and visibility and considered missing
        ageThreshold = 10; %RY need to adjust
        visibilityThreshold = 0.6; %RY need to adjust
        lostThreshold = 10; %RY need to adjust

        % Find the indices of objTracks to delete
        newInds = objTracks.TotalFrameCount <= ageThreshold;
        lowVisibilityInds = visibility < visibilityThreshold;
        lostInds = objTracks.ConsecutiveInvisibleCount >= lostThreshold;
        
        deleteInds = (newInds & lowVisibilityInds) | lostInds;

        % Delete tracks
        objTracks = objTracks(~deleteInds,:);
    end

    % Create new objTracks for unassigned objDetections
    for idx = 1:height(objDetections)
        if ~objDetections.assigned(idx)

            % Set Kalman Filter parameters
            InitialLocation = objDetections.objCentroid(idx,:);
            FilterType = "ConstantVelocity";
            InitialEstimateError = [200,50];
            MotionNoise = [100, 25]; 
            MeasurementNoise = 100;

            % Create a Kalman Filter object
            KalmanFilter = configureKalmanFilter(FilterType,InitialLocation, ...
                InitialEstimateError, MotionNoise, MeasurementNoise);
            
            % Set initial variables for new objTrack
            DetectedLocation = InitialLocation;
            TrackedLocation = InitialLocation;
            PredictedLocation = InitialLocation;
            PredictedState = zeros(1,4);
            TotalFrameCount = 1;
            TotalVisibleCount = 1;
            ConsecutiveInvisibleCount = 0;
            Visible = true;
            Confirmed = false;

            % Convert Kalman Filter to a cell for stacking
            KalmanFilter = {KalmanFilter};

            % Create new objTrack row
            newobjTrack = table(objTrackID,...
                DetectedLocation,...
                PredictedLocation,...
                TrackedLocation,...
                PredictedState,...
                KalmanFilter,...
                TotalFrameCount,...
                TotalVisibleCount,...
                ConsecutiveInvisibleCount,...
                Visible,...
                Confirmed);
            
            % Add it to the array of tracks
            objTracks = [objTracks; newobjTrack];

            % Iterate objTrackID
            objTrackID = objTrackID + 1;
        end
    end
end

%% FUNCTIONS - cornerDetections = detectCorners(pgx,pgy);
function cornerDetections = detectCorners(pgx,pgy)
    minParaPoints = [pgx(1:1:end-1),pgy(1:1:end-1)];
    cornerDetections = table(minParaPoints);
end

%% FUNCTIONS - cornerTracks = predictcornerTracks(cornerTracks);
function cornerTracks = predictcornerTracks(cornerTracks)
    for idx = 1:height(cornerTracks)
        cornerTracks.PredictedLocation(idx,:) = predict(cornerTracks.KalmanFilter{idx});
    end
end

%% FUNCTIONS - [cornerTracks, cornerDetections] = assigncornerTracksDetections(objTracks,cornerTracks,cornerDetections);
function [cornerTracks, cornerDetections] = assigncornerTracksDetections(cornerTracks,cornerDetections)
    % Compute the cost of assigning each corner detection to each corner track
    cost = zeros(height(cornerTracks),height(cornerDetections));
    for idx = 1:height(cornerTracks)
        cost(idx,:) = distance(cornerTracks.KalmanFilter{idx},cornerDetections.minParaPoints);
    end
    
    % Solve the assignment problem
    costOfNonAssignment = 1000;
    assignedIdxPairs = assignDetectionsToTracks(cost,costOfNonAssignment);

    % Set cornerdetected flag & add cornerdetected location to cornerdetected tracks
    if ~isempty(cornerTracks)
        cornerTracks.Visible(:) = false;
        cornerTracks.Visible(assignedIdxPairs(:,1)) = true;
        cornerTracks.DetectedLocation(assignedIdxPairs(:,1),:) = cornerDetections.minParaPoints(assignedIdxPairs(:,2),:);
    end

    % Set assignment flag to corner detections
    if ~isempty(cornerDetections)
        cornerDetections.assigned(:) = false;
        cornerDetections.assigned(assignedIdxPairs(:,2)) = true;
    end
end

%% FUNCTIONS - cornerTracks = updatecornerTracks(cornerTracks,cornerDetections);
function cornerTracks = updatecornerTracks(cornerTracks,cornerDetections)
    
    % Initialize persistent corner counter
    persistent cornerTrackID;
    if isempty(cornerTrackID)
        cornerTrackID = 1;
    end

    % Update tracked corner locations using available corner detections
    for idx = 1:height(cornerTracks)
        if cornerTracks.Visible(idx)
            cornerTracks.TrackedLocation(idx,:) = correct(cornerTracks.KalmanFilter{idx},cornerTracks.DetectedLocation(idx,:));
        end
    end

    % Update remaining corner track variables
    if ~isempty(cornerTracks)
        % Set cornerTracked location to prediction for undetected cornerTracks
        cornerTracks.TrackedLocation(~cornerTracks.Visible,:) = cornerTracks.PredictedLocation(~cornerTracks.Visible,:);

        % Increment all track ages
        cornerTracks.TotalFrameCount = cornerTracks.TotalFrameCount + 1;
        
        % Increment all visible tracks total visible count
        cornerTracks.TotalVisibleCount(cornerTracks.Visible) = cornerTracks.TotalVisibleCount(cornerTracks.Visible) + 1;
        
        % Confirm tracks detected enough
        cornerTrackConfirmationThreshold = 3;
        cornerTracks.Confirmed = cornerTracks.TotalVisibleCount > cornerTrackConfirmationThreshold;

        % Reset the counter for consecutive frames undetected if detected
        cornerTracks.ConsecutiveInvisibleCount(cornerTracks.Visible) = 0;

        % Increment the counter for consecutive frames undetected
        cornerTracks.ConsecutiveInvisibleCount(~cornerTracks.Visible) = cornerTracks.ConsecutiveInvisibleCount(~cornerTracks.Visible) + 1;

        % Compute the fraction of each cornerTracks's age for which it was vis
        visibility = cornerTracks.TotalVisibleCount ./ cornerTracks.TotalFrameCount;
        
        % Set thresholds for age and visibility and considered missing
        ageThreshold = 10; %RY need to adjust
        visibilityThreshold = 0.6; %RY need to adjust
        lostThreshold = 10; %RY need to adjust

        % Find the indices of cornerTracks to delete
        newInds = cornerTracks.TotalFrameCount <= ageThreshold;
        lowVisibilityInds = visibility < visibilityThreshold;
        lostInds = cornerTracks.ConsecutiveInvisibleCount >= lostThreshold;
        
        deleteInds = (newInds & lowVisibilityInds) | lostInds;

        % Delete tracks
        cornerTracks = cornerTracks(~deleteInds,:);
    end

    % Create new cornerTracks for unassigned cornerDetections
    for idx = 1:height(cornerDetections)
        if ~cornerDetections.assigned(idx)

            % Set Kalman Filter parameters
            InitialLocation = cornerDetections.minParaPoints(idx,:);
            FilterType = "ConstantVelocity";
            InitialEstimateError = [200,50];
            MotionNoise = [100, 25]; 
            MeasurementNoise = 100;

            % Create a Kalman Filter object
            KalmanFilter = configureKalmanFilter(FilterType,InitialLocation, ...
                InitialEstimateError, MotionNoise, MeasurementNoise);
            
            % Set initial variables for new cornerTrack
            DetectedLocation = InitialLocation;
            TrackedLocation = InitialLocation;
            PredictedLocation = InitialLocation;
            PredictedState = zeros(1,4);
            TotalFrameCount = 1;
            TotalVisibleCount = 1;
            ConsecutiveInvisibleCount = 0;
            Visible = true;
            Confirmed = false;

            % Convert Kalman Filter to a cell for stacking
            KalmanFilter = {KalmanFilter};

            % Create new cornerTrack row
            newcornerTrack = table(cornerTrackID,...
                DetectedLocation,...
                PredictedLocation,...
                TrackedLocation,...
                PredictedState,...
                KalmanFilter,...
                TotalFrameCount,...
                TotalVisibleCount,...
                ConsecutiveInvisibleCount,...
                Visible,...
                Confirmed);
            
            % Add it to the array of tracks
            cornerTracks = [cornerTracks; newcornerTrack];

            % Iterate objTrackID
            cornerTrackID = cornerTrackID + 1;
        end
    end
end

end
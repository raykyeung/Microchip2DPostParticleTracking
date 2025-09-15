function [myFigures,time_ms,position_x_CS,velocity_x_CS,acceleration_x_CS,position_y_CS,velocity_y_CS,acceleration_y_CS,meanThroughput_x_CS] = B2KParticleTracking30(i,vidObj,myFigures,onePTrackOutputFolderPath,L_start,R_end,tFrames,nFrames,I,CS,x_pos_caption,y_pos_caption,EXP_actualframerate,flagDebug,flagRunWithDisplay)
%% B2KParticleTracking30_pred.m - [Function] Batch single p-Chip tracking
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

%%
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

%=========================================================================%
% function [myFigures, position_x_CS, position_y_CS] = B2KParticleTracking10_pred(i, vidObj, myFigures, onePTrackOutputFolderPath, L_start, R_end, tFrames, nFrames, I, CS, x_pos_caption, y_pos_caption, EXP_actualframerate, flagDebug, flagRunWithDisplay)
% ParticleTracking08_pred processes a video (via vidObj) to track a single
% particle and produces four output videos:
%   - trackedVideo1.mp4: Original (raw) video.
%   - trackedVideo2.mp4: Current-frame detection overlay (yellow mask + blue detection marker).
%   - trackedVideo3.mp4: Accumulated detections overlay (blue markers with low opacity).
%   - trackedVideo5.mp4: Tracking overlay using a constant velocity Kalman filter.
%
% Finally, these videos are concatenated vertically into:
%   B2K_p-Chip_Tracking.mp4
%
% OUTPUTS:
%   myFigures                - Updated figure handles.
%   position_x_CS, position_y_CS - Tracked centroid positions (from constant velocity tracking).

%% Tracking setup
L_start_orig = L_start;
R_end_orig = R_end;

%% Preallocate for position, velocity, acceleration data for each frames
position_x_CS_DETECTED = NaN([1,tFrames],'double');
velocity_x_CS_DETECTED = NaN([1,tFrames],'double');
acceleration_x_CS_DETECTED = NaN([1,tFrames],'double');

position_y_CS_DETECTED = NaN([1,tFrames],'double');
velocity_y_CS_DETECTED = NaN([1,tFrames],'double');
acceleration_y_CS_DETECTED = NaN([1,tFrames],'double');

velocity_y_CS_avg_DETECTED = NaN([1,tFrames],'double');
acceleration_y_CS_avg_DETECTED = NaN([1,tFrames],'double');

position_x_CS_TRACKED = NaN([1,tFrames],'double');
velocity_x_CS_TRACKED = NaN([1,tFrames],'double');
acceleration_x_CS_TRACKED = NaN([1,tFrames],'double');

position_y_CS_TRACKED = NaN([1,tFrames],'double');
velocity_y_CS_TRACKED = NaN([1,tFrames],'double');
acceleration_y_CS_TRACKED = NaN([1,tFrames],'double');

velocity_y_CS_avg_TRACKED = NaN([1,tFrames],'double');
acceleration_y_CS_avg_TRACKED = NaN([1,tFrames],'double');

frame_array = 1:1:tFrames;
frame_num_array = L_start_orig:R_end_orig;
time_ms_array = frame_array*(1/EXP_actualframerate)*(1000);

%% Preallocate blobCentroid and Harris Centroid data (comparing both methods for centroid detection)
blobCentroid_x = NaN([1,tFrames],'double');
blobCentroid_y = NaN([1,tFrames],'double');
HarrisCentroid_x = NaN([1,tFrames],'double');
HarrisCentroid_y = NaN([1,tFrames],'double');

blobCentroid_x_CS_DETECTED = NaN([1,tFrames],'double');
blobCentroid_y_CS_DETECTED = NaN([1,tFrames],'double');
HarrisCentroid_x_CS = NaN([1,tFrames],'double');
HarrisCentroid_y_CS = NaN([1,tFrames],'double');

centroid_distance_CS = NaN([1,tFrames],'double');

%% Video Writer Setup (naming similar to sample code)
vidNumStr = sprintf('_vid%d', i);

% trackedVideo1: Original video output
vidOutputOrigFileName = 'trackedVideo1.mp4';
vidOutputFileName = insertBefore(vidOutputOrigFileName, '.mp4', vidNumStr);
trackedVid1File = fullfile(onePTrackOutputFolderPath, vidOutputFileName);
vw1 = VideoWriter(trackedVid1File, 'MPEG-4');
vw1.FrameRate = 15;
open(vw1);

% trackedVideo2: Current-frame detection overlay
vidOutputOrigFileName = 'trackedVideo2.mp4';
vidOutputFileName = insertBefore(vidOutputOrigFileName, '.mp4', vidNumStr);
trackedVid2File = fullfile(onePTrackOutputFolderPath, vidOutputFileName);
vw2 = VideoWriter(trackedVid2File, 'MPEG-4');
vw2.FrameRate = 15;
open(vw2);

% trackedVideo3: Accumulated detections overlay
vidOutputOrigFileName = 'trackedVideo3.mp4';
vidOutputFileName = insertBefore(vidOutputOrigFileName, '.mp4', vidNumStr);
trackedVid3File = fullfile(onePTrackOutputFolderPath, vidOutputFileName);
vw3 = VideoWriter(trackedVid3File, 'MPEG-4');
vw3.FrameRate = 15;
open(vw3);

% trackedVideo5: Constant Velocity Tracking overlay
vidOutputOrigFileName = 'trackedVideo5.mp4';
vidOutputFileName = insertBefore(vidOutputOrigFileName, '.mp4', vidNumStr);
trackedVid5File = fullfile(onePTrackOutputFolderPath, vidOutputFileName);
vw5 = VideoWriter(trackedVid5File, 'MPEG-4');
vw5.FrameRate = 15;
open(vw5);

%% Preallocation and Setup
position_x_CS_DETECTED = [];  % Will be set after processing
position_y_CS_DETECTED = [];
trackingDataCV = []; % Each row will be [x y] (wrt to intrinsic cooordinate system, px) for frames with tracking data

% Structures for accumulating image data and detections.
utilities.accumulatedImage = [];
utilities.accumulatedDetections = [];

% For constant velocity tracking (CV) overlays.
accumulatedTrackingsCV = [];

% Default parameters (detection and thresholds)
param = getDefaultParameters();
param.entranceThreshold = 20;   % Begin tracking when x > 20. %(x,y) image size is 76 px, 1280 px
param.exitThreshold     = 1260;  % Stop tracking when x > 1260.
param.L_start = L_start;
param.R_end   = R_end;

% CV tracking parameters: use Constant Velocity model.
cvParam = getDefaultParameters();
cvParam.motionModel = 'ConstantVelocity';
cvParam.initialEstimateError = cvParam.initialEstimateError(1:2);
cvParam.motionNoise = cvParam.motionNoise(1:2);
kalmanFilterCV = [];

% Create foreground detector and blob analyzer once.
fgDetector = vision.ForegroundDetector('NumTrainingFrames', 10, 'InitialVariance', param.segmentationThreshold);
blobAnalyzer = vision.BlobAnalysis('MinimumBlobArea',50,'MaximumBlobArea',500,'PerimeterOutputPort',true,'OrientationOutputPort',true,'ExcludeBorderBlobs',true);

% Create a structuring element for mask cleaning.
se = strel('disk', 3);

% Create/update figure handles.
myFigures = B2KMyFigures(myFigures, flagRunWithDisplay);

%% Process Frames (from L_start to R_end)
frameCount = 0;

for frameIdx = L_start:R_end
    frameCount = frameCount + 1;
    frame = read(vidObj, frameIdx);
    frame = ensureRGB(frame);
    frame = im2uint8(frame);  % Ensure frame type for VideoWriter

    %% (1) Write Original Video (trackedVideo1)

    text0originLeftTop = [0 60];
    text1originLeftTop = [0 0];
    text2spacing = [size(frame,1) 0]; %[76 0]
    text2addspacing = [size(frame,1) 0]; %[76 0]
    text3originCenterTop = [(size(frame,2)/2) 0]; %[640 0]

    frameOrig = frame;
    Itext0_DETECTED = insertText(frameOrig,text0originLeftTop,['Captured FPS: ' num2str(EXP_actualframerate) '   |   Playback speed: ' num2str(EXP_actualframerate/vw1.FrameRate,'%.2f') 'X slower'],'AnchorPoint','LeftTop','FontSize',10,'BoxColor','red','BoxOpacity',0.1,'TextColor','white'); %bottom left red text box describing capture FPS & playback in exported video
    Itext1_DETECTED = insertText(Itext0_DETECTED,text1originLeftTop,['Frame #: ' num2str(frameIdx)],'AnchorPoint','LeftTop','FontSize',10,'BoxColor','green','BoxOpacity',0.1,'TextColor','white'); %top left green text box describing Frame # in exported video
    Itext3_DETECTED = insertText(Itext1_DETECTED,text3originCenterTop,'Video with reduced playback speed','AnchorPoint','CenterTop','FontSize',10,'BoxColor','black','BoxOpacity',0.1,'TextColor','white');
    writeVideo(vw1, Itext3_DETECTED);

    %% (2) Obtain Detection and Foreground Mask
    [objArea,objCentroid,bboxOut,isDetected,mask] = detectObjExtended(frame, fgDetector, blobAnalyzer);
    
    numObj = length(objArea); %blob number

    % Update accumulated image and detections (for overlays).
    utilities.accumulatedImage = updateAccumulatedImage(utilities.accumulatedImage, frame);
    if isDetected
        utilities.accumulatedDetections = updateAccumulator(utilities.accumulatedDetections, objCentroid);
    end

    %% (3) Write Current-frame Detection Overlay (trackedVideo2)
    text0originLeftTop = [0 60];
    text1originLeftTop = [0 0];
    text2spacing = [size(frame,1) 0]; %[76 0]
    text2addspacing = [size(frame,1) 0]; %[76 0]
    text3originCenterTop = [(size(frame,2)/2) 0]; %[640 0]

    maskClean = imopen(mask, se);
    frameDet = imoverlay(frame, maskClean > 0, [1, 1, 0]); % yellow overlay for blob pixels
    if isDetected
        frameDet = insertShape(frameDet, 'FilledCircle', [objCentroid, 3], 'Color', 'blue', 'Opacity', 0.6);
    end

    Itext0_DETECTED = insertText(frameDet,text0originLeftTop,['Captured FPS: ' num2str(EXP_actualframerate) '   |   Playback speed: ' num2str(EXP_actualframerate/vw1.FrameRate,'%.2f') 'X slower'],'AnchorPoint','LeftTop','FontSize',10,'BoxColor','red','BoxOpacity',0.1,'TextColor','white'); %bottom left red text box describing capture FPS & playback in exported video
    Itext1_DETECTED = insertText(Itext0_DETECTED,text1originLeftTop,['Frame #: ' num2str(frameIdx)],'AnchorPoint','LeftTop','FontSize',10,'BoxColor','green','BoxOpacity',0.1,'TextColor','white'); %top left green text box describing Frame # in exported video
    
    if ~isDetected
        Itext2_DETECTED = insertText(Itext1_DETECTED,text2spacing,'No p-Chip found','AnchorPoint','LeftTop','FontSize',10,'BoxColor','blue','BoxOpacity',0.2,'TextColor','white'); %top left blue text box describing p-Chip detection/tracking information in exported video
    else
        % Transform data to proper coordinate system
        [objCentroid_x_CS,objCentroid_y_CS] = intrinsicToWorld(CS,objCentroid(1,1),objCentroid(1,2));

        % Store objCentroid data in blobCentroid array
        blobCentroid_x_CS_DETECTED(frameIdx - L_start_orig + 1) = objCentroid_x_CS;
        blobCentroid_y_CS_DETECTED(frameIdx - L_start_orig + 1) = objCentroid_y_CS;

        % Position calculation [mm]
        position_x_CS_DETECTED(frameIdx - L_start_orig + 1) = objCentroid_x_CS;
        position_y_CS_DETECTED(frameIdx - L_start_orig + 1) = objCentroid_y_CS;

        firstframe_DETECTED = frame_num_array(find(~isnan(position_x_CS_DETECTED),1));

        if frameIdx == firstframe_DETECTED
            velocity_x_CS_DETECTED(frameIdx - L_start_orig + 1) = NaN;
            acceleration_x_CS_DETECTED(frameIdx - L_start_orig + 1) = NaN;
            velocity_y_CS_DETECTED(frameIdx - L_start_orig + 1) = NaN;
            acceleration_y_CS_DETECTED(frameIdx - L_start_orig + 1) = NaN;
        elseif frameIdx == firstframe_DETECTED + 1
            velocity_x_CS_DETECTED(frameIdx - L_start_orig + 1) = (position_x_CS_DETECTED(frameIdx - L_start_orig + 1) - position_x_CS_DETECTED(frameIdx - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]
            acceleration_x_CS_DETECTED(frameIdx - L_start_orig + 1) = NaN;
            velocity_y_CS_DETECTED(frameIdx - L_start_orig + 1) = (position_y_CS_DETECTED(frameIdx - L_start_orig + 1) - position_y_CS_DETECTED(frameIdx - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]
            acceleration_y_CS_DETECTED(frameIdx - L_start_orig + 1) = NaN;
        else
            % Velocity calculation [m/s]
            velocity_x_CS_DETECTED(frameIdx - L_start_orig + 1) = (position_x_CS_DETECTED(frameIdx - L_start_orig + 1) - position_x_CS_DETECTED(frameIdx - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]
            velocity_y_CS_DETECTED(frameIdx - L_start_orig + 1) = (position_y_CS_DETECTED(frameIdx - L_start_orig + 1) - position_y_CS_DETECTED(frameIdx - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]

            % Acceleration calculation [m/s^2]
            acceleration_x_CS_DETECTED(frameIdx - L_start_orig + 1) = (velocity_x_CS_DETECTED(frameIdx - L_start_orig + 1) - velocity_x_CS_DETECTED(frameIdx - L_start_orig))*(EXP_actualframerate); %[m/s^2]
            acceleration_y_CS_DETECTED(frameIdx - L_start_orig + 1) = (velocity_y_CS_DETECTED(frameIdx - L_start_orig + 1) - velocity_y_CS_DETECTED(frameIdx - L_start_orig))*(EXP_actualframerate); %[m/s^2]
        end

        %Note: Take movmean for velocity either after the end or do the movmean at real time 
        % Multi-point average of velocity on x & y by 5 points
        velocity_x_CS_avg_DETECTED = movmean(velocity_x_CS_DETECTED, 5);
        velocity_y_CS_avg_DETECTED = movmean(velocity_y_CS_DETECTED, 5);

        % Multi-point average of acceleration on x & y by 5 points
        acceleration_x_CS_avg_DETECTED = movmean(acceleration_x_CS_DETECTED, 5);
        acceleration_y_CS_avg_DETECTED = movmean(acceleration_y_CS_DETECTED, 5);

        text_str = cell(numObj,1);
        text2origin = zeros(numObj,2);
    
        % [Old iterative figure plotting]
        for ii = 1:numObj
            text_str{ii} = ['[Detected p-Chip' num2str(ii) ']' x_pos_caption num2str(objCentroid_x_CS,'%.2f') ',' y_pos_caption num2str(objCentroid_y_CS,'%.2f')];
            text2origin(ii,:) = text1originLeftTop + text2spacing;
            text2spacing = text2spacing + text2addspacing;
        end
        Itext2_DETECTED = insertText(Itext1_DETECTED,text2origin,text_str,'AnchorPoint','LeftTop','FontSize',10,'BoxColor','blue','BoxOpacity',0.2,'TextColor','white');
    end

    %Display image
    if flagRunWithDisplay == 1
        imshow(Itext2_DETECTED,'InitialMagnification',500); %Can suppress this to improve performance
    end

    Itext3_DETECTED = insertText(Itext2_DETECTED,text3originCenterTop,'Gaussian Mixture Model with Blob Analysis (blob: yellow, blob centroid: blue circle)','AnchorPoint','CenterTop','FontSize',10,'BoxColor','black','BoxOpacity',0.1,'TextColor','white');
    writeVideo(vw2, Itext3_DETECTED);

    %% (4) Write Accumulated Detections Overlay (trackedVideo3)
    text0originLeftTop = [0 60];
    text1originLeftTop = [0 0];
    text2spacing = [size(frame,1) 0]; %[76 0]
    text2addspacing = [size(frame,1) 0]; %[76 0]
    text3originCenterTop = [(size(frame,2)/2) 0]; %[640 0]

    frameAccum = frame;
    if ~isempty(utilities.accumulatedDetections)
        frameAccum = insertShape(frameAccum, 'FilledCircle', [utilities.accumulatedDetections, repmat(3, size(utilities.accumulatedDetections,1), 1)], 'Color', 'blue', 'Opacity', 0.05);
    end

    Itext0_DETECTED = insertText(frameAccum,text0originLeftTop,['Captured FPS: ' num2str(EXP_actualframerate) '   |   Playback speed: ' num2str(EXP_actualframerate/vw1.FrameRate,'%.2f') 'X slower'],'AnchorPoint','LeftTop','FontSize',10,'BoxColor','red','BoxOpacity',0.1,'TextColor','white'); %bottom left red text box describing capture FPS & playback in exported video
    Itext1_DETECTED = insertText(Itext0_DETECTED,text1originLeftTop,['Frame #: ' num2str(frameIdx)],'AnchorPoint','LeftTop','FontSize',10,'BoxColor','green','BoxOpacity',0.1,'TextColor','white'); %top left green text box describing Frame # in exported video
 
    if ~isDetected
        Itext2_DETECTED = insertText(Itext1_DETECTED,text2spacing,'No p-Chip found','AnchorPoint','LeftTop','FontSize',10,'BoxColor','blue','BoxOpacity',0.2,'TextColor','white'); %top left blue text box describing p-Chip detection/tracking information in exported video
    else
        text_str = cell(numObj,1);
        text2origin = zeros(numObj,2);
    
        % [Old iterative figure plotting]
        for ii = 1:numObj
            text_str{ii} = ['[Detected p-Chip' num2str(ii) ']' x_pos_caption num2str(objCentroid_x_CS,'%.2f') ',' y_pos_caption num2str(objCentroid_y_CS,'%.2f')];
            text2origin(ii,:) = text1originLeftTop + text2spacing;
            text2spacing = text2spacing + text2addspacing;
        end
        Itext2_DETECTED = insertText(Itext1_DETECTED,text2origin,text_str,'AnchorPoint','LeftTop','FontSize',10,'BoxColor','blue','BoxOpacity',0.2,'TextColor','white');    
    end
    Itext3_DETECTED = insertText(Itext2_DETECTED,text3originCenterTop,'Accumulated detections (blue circles)','AnchorPoint','CenterTop','FontSize',10,'BoxColor','black','BoxOpacity',0.1,'TextColor','white');
    
    % frameAccum = insertText(frameAccum, [10, 10], 'Accumulated Detections: Blue (Opacity 0.05)', 'FontSize', 12, 'BoxColor', 'black', 'TextColor', 'white', 'BoxOpacity', 0.6);
    writeVideo(vw3, Itext3_DETECTED);

    %% (5) Process Constant Velocity Tracking for trackedVideo5
    text0originLeftTop = [0 60];
    text1originLeftTop = [0 0];
    text2spacing = [size(frame,1) 0]; %[76 0]
    text2addspacing = [size(frame,1) 0]; %[76 0]
    text3originCenterTop = [(size(frame,2)/2) 0]; %[640 0]

    trackedCV = [];
    labelCV = '';
    if isDetected && objCentroid(1) >= param.entranceThreshold
        if isempty(kalmanFilterCV)
            initLoc = computeInitialLocation(cvParam, objCentroid);
            kalmanFilterCV = configureKalmanFilter(cvParam.motionModel, initLoc, cvParam.initialEstimateError, cvParam.motionNoise, cvParam.measurementNoise);
            trackedCV = correct(kalmanFilterCV, objCentroid);
            [trackedCV_x_CS, trackedCV_y_CS] = intrinsicToWorld(CS,trackedCV(1,1),trackedCV(1,2));
            trackedCV_CS = [trackedCV_x_CS, trackedCV_y_CS];
            labelCV = 'Initial';
        else
            predict(kalmanFilterCV);
            trackedCV = correct(kalmanFilterCV, objCentroid);
            [trackedCV_x_CS, trackedCV_y_CS] = intrinsicToWorld(CS,trackedCV(1,1),trackedCV(1,2));
            trackedCV_CS = [trackedCV_x_CS, trackedCV_y_CS];
            labelCV = 'Corrected';
        end
    end
    if ~isempty(trackedCV) && trackedCV(1) > param.exitThreshold
        trackedCV = [];
        labelCV = '';
    end
    if ~isempty(trackedCV)
        accumulatedTrackingsCV = updateAccumulator(accumulatedTrackingsCV, trackedCV);
        trackingDataCV = [trackingDataCV; trackedCV];
    end

    % Build overlay for constant velocity tracking (trackedVideo5)
    overlayCV = frame;
    % Overlay accumulated detections.
    if ~isempty(utilities.accumulatedDetections)
        overlayCV = insertShape(overlayCV, 'FilledCircle', [utilities.accumulatedDetections, repmat(3, size(utilities.accumulatedDetections,1), 1)], 'Color', 'blue', 'Opacity', 0.6);
    end
    % Overlay accumulated tracking markers.
    if ~isempty(accumulatedTrackingsCV)
        overlayCV = insertShape(overlayCV, 'Circle', [accumulatedTrackingsCV, repmat(4, size(accumulatedTrackingsCV,1), 1)], 'Color', 'red', 'LineWidth', 1);
        if size(accumulatedTrackingsCV,1) >= 2
            segs = [accumulatedTrackingsCV(1:end-1,:) accumulatedTrackingsCV(2:end,:)];
            overlayCV = insertShape(overlayCV, 'Line', segs, 'Color', 'red', 'LineWidth', 1);
        end
    end
    % Overlay current detection and tracking markers.
    if isDetected
        overlayCV = insertShape(overlayCV, 'FilledCircle', [objCentroid, 3], 'Color', 'blue', 'Opacity', 0.6);
    end
    if ~isempty(trackedCV)
        overlayCV = insertShape(overlayCV, 'Circle', [trackedCV, 4], 'Color', 'red', 'LineWidth', 1);
        overlayCV = insertText(overlayCV, trackedCV + [10, 10], labelCV, 'FontSize', 12, 'BoxColor', 'red', 'TextColor', 'white', 'BoxOpacity', 0.1);
    end
    % legendText = 'Detection: blue filled circle, Tracking: red circle/line';
    % overlayCV = insertText(overlayCV, [10, 10], legendText, 'FontSize', 12, 'BoxColor', 'black', 'TextColor', 'white', 'BoxOpacity', 0.1);

    Itext0_DETECTED = insertText(overlayCV,text0originLeftTop,['Captured FPS: ' num2str(EXP_actualframerate) '   |   Playback speed: ' num2str(EXP_actualframerate/vw1.FrameRate,'%.2f') 'X slower'],'AnchorPoint','LeftTop','FontSize',10,'BoxColor','red','BoxOpacity',0.1,'TextColor','white'); %bottom left red text box describing capture FPS & playback in exported video
    Itext1_DETECTED = insertText(Itext0_DETECTED,text1originLeftTop,['Frame #: ' num2str(frameIdx)],'AnchorPoint','LeftTop','FontSize',10,'BoxColor','green','BoxOpacity',0.1,'TextColor','white'); %top left green text box describing Frame # in exported video
 
    if ~isDetected && isempty(trackedCV)
        Itext2_DETECTED = insertText(Itext1_DETECTED,text2spacing,'p-Chip not detected or tracked','AnchorPoint','LeftTop','FontSize',10,'BoxColor','blue','BoxOpacity',0.2,'TextColor','white'); %top left blue text box describing p-Chip detection/tracking information in exported video
    elseif isDetected && isempty(trackedCV)
        text_str = cell(numObj,1);
        text2origin = zeros(numObj,2);
    
        % [Old iterative figure plotting]
        for ii = 1:numObj
            text_str{ii} = ['[Detected p-Chip' num2str(ii) ']' x_pos_caption num2str(objCentroid_x_CS,'%.2f') ',' y_pos_caption num2str(objCentroid_y_CS,'%.2f')];
            text2origin(ii,:) = text1originLeftTop + text2spacing;
            text2spacing = text2spacing + text2addspacing;
        end
        Itext2_DETECTED = insertText(Itext1_DETECTED,text2origin,text_str,'AnchorPoint','LeftTop','FontSize',10,'BoxColor','blue','BoxOpacity',0.2,'TextColor','white');    
    else %isDetected && ~isempty(trackedCV)
        text_str = cell(numObj,1);
        text2origin = zeros(numObj,2);
    
        % [Old iterative figure plotting]
        for ii = 1:numObj
            text_str{ii} = ['[Detected and tracked p-Chip' num2str(ii) ']' x_pos_caption num2str(trackedCV_x_CS,'%.2f') ',' y_pos_caption num2str(trackedCV_y_CS,'%.2f')];
            text2origin(ii,:) = text1originLeftTop + text2spacing;
            text2spacing = text2spacing + text2addspacing;
        end
        Itext2_DETECTED = insertText(Itext1_DETECTED,text2origin,text_str,'AnchorPoint','LeftTop','FontSize',10,'BoxColor','blue','BoxOpacity',0.2,'TextColor','white');    
    end
    Itext3_DETECTED = insertText(Itext2_DETECTED,text3originCenterTop,'p-Chip tracking with Kalman filter prediction (detections: blue, tracks: red)','AnchorPoint','CenterTop','FontSize',10,'BoxColor','black','BoxOpacity',0.1,'TextColor','white');

    writeVideo(vw5, Itext3_DETECTED);

    %%
    if isempty(trackedCV)

    else
        % Transform data to proper coordinate system
        [trackedCV_x_CS,trackedCV_y_CS] = intrinsicToWorld(CS,trackedCV(1,1),trackedCV(1,2));

        % Store objCentroid data in blobCentroid array
        blobCentroid_x_CS_TRACKED(frameIdx - L_start_orig + 1) = trackedCV_x_CS;
        blobCentroid_y_CS_TRACKED(frameIdx - L_start_orig + 1) = trackedCV_y_CS;

        % Position calculation [mm]
        position_x_CS_TRACKED(frameIdx - L_start_orig + 1) = trackedCV_x_CS;
        position_y_CS_TRACKED(frameIdx - L_start_orig + 1) = trackedCV_y_CS;

        firstframe_TRACKED = frame_num_array(find(~isnan(position_x_CS_TRACKED),1));

        if frameIdx == firstframe_TRACKED
            velocity_x_CS_TRACKED(frameIdx - L_start_orig + 1) = NaN;
            acceleration_x_CS_TRACKED(frameIdx - L_start_orig + 1) = NaN;
            velocity_y_CS_TRACKED(frameIdx - L_start_orig + 1) = NaN;
            acceleration_y_CS_TRACKED(frameIdx - L_start_orig + 1) = NaN;
        elseif frameIdx == firstframe_TRACKED + 1
            velocity_x_CS_TRACKED(frameIdx - L_start_orig + 1) = (position_x_CS_TRACKED(frameIdx - L_start_orig + 1) - position_x_CS_TRACKED(frameIdx - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]
            acceleration_x_CS_TRACKED(frameIdx - L_start_orig + 1) = NaN;
            velocity_y_CS_TRACKED(frameIdx - L_start_orig + 1) = (position_y_CS_TRACKED(frameIdx - L_start_orig + 1) - position_y_CS_TRACKED(frameIdx - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]
            acceleration_y_CS_TRACKED(frameIdx - L_start_orig + 1) = NaN;
        else
            % Velocity calculation [m/s]
            velocity_x_CS_TRACKED(frameIdx - L_start_orig + 1) = (position_x_CS_TRACKED(frameIdx - L_start_orig + 1) - position_x_CS_TRACKED(frameIdx - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]
            velocity_y_CS_TRACKED(frameIdx - L_start_orig + 1) = (position_y_CS_TRACKED(frameIdx - L_start_orig + 1) - position_y_CS_TRACKED(frameIdx - L_start_orig))*(EXP_actualframerate)/1000; %[m/s]

            % Acceleration calculation [m/s^2]
            acceleration_x_CS_TRACKED(frameIdx - L_start_orig + 1) = (velocity_x_CS_TRACKED(frameIdx - L_start_orig + 1) - velocity_x_CS_TRACKED(frameIdx - L_start_orig))*(EXP_actualframerate); %[m/s^2]
            acceleration_y_CS_TRACKED(frameIdx - L_start_orig + 1) = (velocity_y_CS_TRACKED(frameIdx - L_start_orig + 1) - velocity_y_CS_TRACKED(frameIdx - L_start_orig))*(EXP_actualframerate); %[m/s^2]
        end

        %Note: Take movmean for velocity either after the end or do the movmean at real time 
        % Multi-point average of velocity on x & y by 5 points
        velocity_x_CS_avg_TRACKED = movmean(velocity_x_CS_TRACKED, 5);
        velocity_y_CS_avg_TRACKED = movmean(velocity_y_CS_TRACKED, 5);

        % Multi-point average of acceleration on x & y by 5 points
        acceleration_x_CS_avg_TRACKED = movmean(acceleration_x_CS_TRACKED, 5);
        acceleration_y_CS_avg_TRACKED = movmean(acceleration_y_CS_TRACKED, 5);
    end
end

%% Set Output Tracking Data
if ~isempty(trackingDataCV)
    position_x = trackingDataCV(:, 1);
    position_y = trackingDataCV(:, 2);
else
    position_x = [];
    position_y = [];
end

%% x-position plot
myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);
figure(myFigures.fig.(myFigures.handle(myFigures.figNum)))
figXPosData = gcf;
tileXPosData = tiledlayout(figXPosData,3,1,'TileSpacing','tight','padding','tight');
tileXPosData.TileIndexing = 'rowmajor';
tileXPosData.OuterPosition = [0.25 0 0.5 1];
TitleXPosData = title(tileXPosData,'Position Data','Interpreter','none');
TitleXPosData.FontSize = 30;

%% Remove NaN values at start and end of the array
firstIndex = find(~isnan(position_x_CS_TRACKED), 1, 'first');
lastIndex  = find(~isnan(position_x_CS_TRACKED), 1, 'last');

position_x_CS_TRACKED_trimmed  = position_x_CS_TRACKED(firstIndex:lastIndex);
velocity_x_CS_TRACKED_trimmed = velocity_x_CS_TRACKED(firstIndex:lastIndex);
acceleration_x_CS_TRACKED_trimmed = acceleration_x_CS_TRACKED(firstIndex:lastIndex);

position_y_CS_TRACKED_trimmed  = position_y_CS_TRACKED(firstIndex:lastIndex);
velocity_y_CS_TRACKED_trimmed = velocity_y_CS_TRACKED(firstIndex:lastIndex);
acceleration_y_CS_TRACKED_trimmed = acceleration_y_CS_TRACKED(firstIndex:lastIndex);

%% Correct time_ms_array size to account for reduced analysis window due to entrance and exit threshold
time_ms_array_trimmed = time_ms_array(:,1:size(position_x_CS_TRACKED_trimmed,2));

%% X-Position, -Velocity, -Acceleration Plot
nexttile(1)
sObj_xPos = scatter(time_ms_array_trimmed,position_x_CS_TRACKED_trimmed,'MarkerEdgeColor',[0.1961, 0.5333, 0.7412],'MarkerFaceColor','none','MarkerEdgeAlpha',0.25,'MarkerFaceAlpha',1,'SizeData',18); %SizeData = 36 (default)
txt = title('X-position');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Position [mm]')
hold on

nexttile(1)
position_x_CS_TRACKED_smoothedSG = smoothdata(position_x_CS_TRACKED,'sgolay'); %Savitzy-Golay filter works but provides erroneous smoothing at ends of data (corresponding to NaN values); 'nanflag' name-value argument of 'includenan' does not work
idx_NaN = isnan(position_x_CS_TRACKED);
position_x_CS_TRACKED_smoothedSG(idx_NaN) = NaN; %Make NaN (remove erroneous ends)
pObj_xPos_SG = plot(time_ms_array_trimmed, position_x_CS_TRACKED_smoothedSG(firstIndex:lastIndex), 'Color',[0.8353, 0.2431, 0.3098 1],'Linewidth',2); %[m/s^2]
hold on

nexttile(2)
sObj_xVel = scatter(time_ms_array_trimmed,velocity_x_CS_TRACKED_trimmed,'MarkerEdgeColor',[0.1961, 0.5333, 0.7412],'MarkerFaceColor',[0.1961, 0.5333, 0.7412],'MarkerEdgeAlpha',0.25,'MarkerFaceAlpha',0.1,'SizeData',18);
txt = title('X-velocity');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Velocity [m/s]')
hold on

nexttile(2)
velocity_x_CS_TRACKED_smoothedSG = smoothdata(velocity_x_CS_TRACKED,'sgolay'); %Savitzy-Golay filter works but provides erroneous smoothing at ends of data (corresponding to NaN values); 'nanflag' name-value argument of 'includenan' does not work
idx_NaN = isnan(velocity_x_CS_TRACKED);

% find the first and last zero
x_firstZero = find(idx_NaN==0, 1, 'first');
x_lastZero = find(idx_NaN==0, 1, 'last');

% flip them to 1
idx_NaN([x_firstZero,x_lastZero]) = 1;

velocity_x_CS_TRACKED_smoothedSG(idx_NaN) = NaN; %Make NaN (remove erroneous ends)
pObj_xVel_SG = plot(time_ms_array_trimmed, velocity_x_CS_TRACKED_smoothedSG(firstIndex:lastIndex), 'Color',[0.8353, 0.2431, 0.3098 1],'Linewidth',2); %[m/s^2]
hold on

nexttile(3)
sObj_xAcc = scatter(time_ms_array_trimmed,acceleration_x_CS_TRACKED_trimmed,'MarkerEdgeColor',[0.1961, 0.5333, 0.7412],'MarkerFaceColor',[0.1961, 0.5333, 0.7412],'MarkerEdgeAlpha',0.25,'MarkerFaceAlpha',0.1,'SizeData',18);
txt = title('X-acceleration');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Acceleration [m/s^{2}]')
hold on

nexttile(3)
acceleration_x_CS_TRACKED_smoothedSG = smoothdata(acceleration_x_CS_TRACKED,'sgolay'); %Savitzy-Golay filter works but provides erroneous smoothing at ends of data (corresponding to NaN values); 'nanflag' name-value argument of 'includenan' does not work
idx_NaN = isnan(acceleration_x_CS_TRACKED);

% find the first and last zero
x_firstZero = find(idx_NaN==0, 2, 'first');
x_lastZero = find(idx_NaN==0, 2, 'last');

% flip them to 1
idx_NaN([x_firstZero,x_lastZero]) = 1;

acceleration_x_CS_TRACKED_smoothedSG(idx_NaN) = NaN; %Make NaN (remove erroneous ends)
pObj_xAcc_SG = plot(time_ms_array_trimmed, acceleration_x_CS_TRACKED_smoothedSG(firstIndex:lastIndex), 'Color',[0.8353, 0.2431, 0.3098 1],'Linewidth',2); %[m/s^2]
hold on

%% Y-Position, -Velocity, -Acceleration Plot

myFigures = B2KMyFigures(myFigures,flagRunWithDisplay);
figure(myFigures.fig.(myFigures.handle(myFigures.figNum)))
figYPosData = gcf;
tileYPosData = tiledlayout(figYPosData,3,1,'TileSpacing','tight','padding','tight');
tileYPosData.TileIndexing = 'rowmajor';
tileYPosData.OuterPosition = [0.25 0 0.5 1];
TitleYPosData = title(tileYPosData,'Position Data','Interpreter','none');
TitleYPosData.FontSize = 30;

nexttile(1)
sObj_yPos = scatter(time_ms_array_trimmed,position_y_CS_TRACKED_trimmed,'MarkerEdgeColor',[0.1961, 0.5333, 0.7412],'MarkerFaceColor','none','MarkerEdgeAlpha',0.25,'MarkerFaceAlpha',1,'SizeData',18); %SizeData = 36 (default)
txt = title('Y-position');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Position [mm]')
hold on

nexttile(1)
position_y_CS_TRACKED_smoothedSG = smoothdata(position_y_CS_TRACKED,'sgolay'); %Savitzy-Golay filter works but provides erroneous smoothing at ends of data (corresponding to NaN values); 'nanflag' name-value argument of 'includenan' does not work
idx_NaN = isnan(position_y_CS_TRACKED);
position_y_CS_TRACKED_smoothedSG(idx_NaN) = NaN; %Make NaN (remove erroneous ends)
pObj_yPos_SG = plot(time_ms_array_trimmed, position_y_CS_TRACKED_smoothedSG(firstIndex:lastIndex), 'Color',[0.8353, 0.2431, 0.3098 1],'Linewidth',2); %[m/s^2]
hold on

nexttile(2)
sObj_yVel = scatter(time_ms_array_trimmed,velocity_y_CS_TRACKED_trimmed,'MarkerEdgeColor',[0.1961, 0.5333, 0.7412],'MarkerFaceColor',[0.1961, 0.5333, 0.7412],'MarkerEdgeAlpha',0.25,'MarkerFaceAlpha',0.1,'SizeData',18);
txt = title('Y-velocity');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Velocity [m/s]')
hold on

nexttile(2)
velocity_y_CS_TRACKED_smoothedSG = smoothdata(velocity_y_CS_TRACKED,'sgolay'); %Savitzy-Golay filter works but provides erroneous smoothing at ends of data (corresponding to NaN values); 'nanflag' name-value argument of 'includenan' does not work
idx_NaN = isnan(velocity_y_CS_TRACKED);

% find the first and last zero
y_firstZero = find(idx_NaN==0, 1, 'first');
y_lastZero = find(idx_NaN==0, 1, 'last');

% flip them to 1
idx_NaN([y_firstZero,y_lastZero]) = 1;

velocity_y_CS_TRACKED_smoothedSG(idx_NaN) = NaN; %Make NaN (remove erroneous ends)
pObj_yVel_SG = plot(time_ms_array_trimmed, velocity_y_CS_TRACKED_smoothedSG(firstIndex:lastIndex), 'Color',[0.8353, 0.2431, 0.3098 1],'Linewidth',2); %[m/s^2]
hold on

nexttile(3)
sObj_yAcc = scatter(time_ms_array_trimmed,acceleration_y_CS_TRACKED_trimmed,'MarkerEdgeColor',[0.1961, 0.5333, 0.7412],'MarkerFaceColor',[0.1961, 0.5333, 0.7412],'MarkerEdgeAlpha',0.25,'MarkerFaceAlpha',0.1,'SizeData',18);
txt = title('Y-acceleration');
txt.FontSize = 16;
txt.Units = "normalized";
set(get(gca,"Title"),"Position",[-0.25 0.5])
xlabel('ms')
ylabel('Acceleration [m/s^{2}]')
hold on

nexttile(3)
acceleration_y_CS_TRACKED_smoothedSG = smoothdata(acceleration_y_CS_TRACKED,'sgolay'); %Savitzy-Golay filter works but provides erroneous smoothing at ends of data (corresponding to NaN values); 'nanflag' name-value argument of 'includenan' does not work
idx_NaN = isnan(acceleration_y_CS_TRACKED);

% find the first and last zero
y_firstZero = find(idx_NaN==0, 2, 'first');
y_lastZero = find(idx_NaN==0, 2, 'last');

% flip them to 1
idx_NaN([y_firstZero,y_lastZero]) = 1;

acceleration_y_CS_TRACKED_smoothedSG(idx_NaN) = NaN; %Make NaN (remove erroneous ends)
pObj_yAcc_SG = plot(time_ms_array_trimmed, acceleration_y_CS_TRACKED_smoothedSG(firstIndex:lastIndex), 'Color',[0.8353, 0.2431, 0.3098 1],'Linewidth',2); %[m/s^2]
hold on

%% Save outputs
time_ms = time_ms_array_trimmed;

position_x_CS = position_x_CS_TRACKED_smoothedSG(firstIndex:lastIndex);
velocity_x_CS = velocity_x_CS_TRACKED_smoothedSG(firstIndex:lastIndex);
acceleration_x_CS = acceleration_x_CS_TRACKED_smoothedSG(firstIndex:lastIndex);

position_y_CS = position_y_CS_TRACKED_smoothedSG(firstIndex:lastIndex);
velocity_y_CS = velocity_y_CS_TRACKED_smoothedSG(firstIndex:lastIndex);
acceleration_y_CS = acceleration_y_CS_TRACKED_smoothedSG(firstIndex:lastIndex);

meanThroughput_x_CS = mean(velocity_x_CS,"omitnan");

%% Close all Video Writers
close(vw1); close(vw2); close(vw3); close(vw5);

%% Concatenate all trackedVideo*.mp4 files into final video.
% finalOutputFile = fullfile(onePTrackOutputFolderPath, 'B2K_p-Chip_Tracking.mp4');
% concatenateTrackedVideos(onePTrackOutputFolderPath, finalOutputFile, 0.25);

%% Display Final Trajectory
myFigures = showTrajectoryFigure(myFigures,flagRunWithDisplay,utilities.accumulatedImage, utilities.accumulatedDetections, accumulatedTrackingsCV);

%%
clear B2KMyFigures

end

%% ===== Local Functions =====

function [objArea,objCentroid,bboxOut,isDetected,mask] = detectObjExtended(frame, fgDetector, blobAnalyzer)
    % detectObjExtended uses the provided foreground detector and blob analyzer
    % to detect the particle in the current frame.
    gray = im2gray(im2single(frame));
    mask = step(fgDetector, gray);
    [objArea,objCentroid,bboxOut,perimeter,orientation] = step(blobAnalyzer, mask);
    if isempty(objCentroid)
        objCentroid = [];
        isDetected = false;
    else
        objCentroid = objCentroid(1, :);  % Use the first detected blob.
        isDetected = true;
    end
end

function accumulator = updateAccumulator(accumulator, newMarker)
    % updateAccumulator appends new marker(s) if available.
    if ~isempty(newMarker)
        accumulator = [accumulator; newMarker];
    end
end

function accumulatedImage = updateAccumulatedImage(accumulatedImage, frame)
    % updateAccumulatedImage computes the pixelwise maximum over frames.
    if isempty(accumulatedImage)
        accumulatedImage = im2single(frame);
    else
        accumulatedImage = max(accumulatedImage, im2single(frame));
    end
end

function concatenateTrackedVideos(folderPath, outputVideoFile, blackFrameFactor)
    % concatenateTrackedVideos finds all files starting with 'trackedVideo' and ending with .mp4
    % in folderPath and concatenates them vertically with a black frame between each pair.
    % The black frame height is specified as a fraction (blackFrameFactor) of the video height.
    %
    % Usage:
    %   concatenateTrackedVideos(folderPath, outputVideoFile)
    %   concatenateTrackedVideos(folderPath, outputVideoFile, blackFrameFactor)
    %
    % Example:
    %   concatenateTrackedVideos('C:\Videos', 'combinedOutput.mp4', 0.25)

    if nargin < 3
        blackFrameFactor = 0.25; % Default black frame height is 25% of the video height
    end

    % Get list of tracked video files
    videoFiles = dir(fullfile(folderPath, 'trackedVideo*.mp4'));
    if isempty(videoFiles)
        error('No tracked video files found in %s.', folderPath);
    end

    % Sort files by name
    [~, idx] = sort({videoFiles.name});
    videoFiles = videoFiles(idx);
    numVideos = numel(videoFiles);

    % Create VideoReader objects for each video file
    vidReaders = cell(1, numVideos);
    for i = 1:numVideos
        vidReaders{i} = VideoReader(fullfile(folderPath, videoFiles(i).name));
    end

    % Read the first frame of the first video to determine frame size
    firstFrame = readFrame(vidReaders{1});
    [frameHeight, frameWidth, numChannels] = size(firstFrame);
    % Reset the first video to the beginning
    vidReaders{1}.CurrentTime = 0;

    % Calculate black frame height in pixels
    blackFrameHeight = round(frameHeight * blackFrameFactor);
    % Create a black frame (all zeros) with the same width and channels as the video frames
    blackFrame = zeros(blackFrameHeight, frameWidth, numChannels, 'uint8');

    % Initialize the output video writer
    outputVideo = VideoWriter(outputVideoFile, 'MPEG-4');
    outputVideo.FrameRate = vidReaders{1}.FrameRate;
    open(outputVideo);

    % Concatenate frames vertically with black spacing in between
    while all(cellfun(@(v) hasFrame(v), vidReaders))
        % Read one frame from each video
        frames = cell(1, numVideos);
        for i = 1:numVideos
            frames{i} = readFrame(vidReaders{i});
        end

        % Start with the first frame then insert a black frame before each subsequent frame
        combinedFrame = frames{1};
        for i = 2:numVideos
            combinedFrame = vertcat(combinedFrame, blackFrame, frames{i});
        end

        % Write the combined frame to the output video
        writeVideo(outputVideo, combinedFrame);
    end

    close(outputVideo);
    fprintf('Combined video written as "%s".\n', outputVideoFile);
end

function rgbFrame = ensureRGB(frame)
    % ensureRGB ensures the frame has three color channels.
    if ndims(frame) == 2
        rgbFrame = repmat(frame, [1, 1, 3]);
    else
        rgbFrame = frame;
    end
end

function out = imoverlay(in, mask, color)
    % imoverlay overlays a binary mask onto an image with the specified color.
    if ~isa(in, 'uint8')
        in = im2uint8(in);
    end
    out = in;
    if isa(in, 'uint8')
        rgbFill = uint8(255 * color);
    else
        rgbFill = color;
    end
    for i = 1:3
        channel = out(:, :, i);
        channel(mask) = rgbFill(i);
        out(:, :, i) = channel;
    end
end

function myFigures = showTrajectoryFigure(myFigures,flagRunWithDisplay,accumulatedImage, accumDets, accumTracks)
    % showTrajectoryFigure displays the accumulated trajectory.
    myFigures = B2KMyFigures(myFigures, flagRunWithDisplay);
    imshow(im2uint8(mat2gray(accumulatedImage))); hold on;
    if ~isempty(accumDets)
        plot(accumDets(:,1), accumDets(:,2), 'b.', 'MarkerSize', 6);
    end
    if ~isempty(accumTracks)
        plot(accumTracks(:,1), accumTracks(:,2), 'ro-', 'LineWidth', 1);
        legend('Detection', 'Tracking');
    end
    title('Accumulated Trajectory');
end

function loc = computeInitialLocation(param, detectedLocation)
    % computeInitialLocation returns the initial location for the Kalman filter.
    if strcmp(param.initialLocation, 'Same as first detection')
        loc = detectedLocation;
    else
        loc = param.initialLocation;
    end
end

function param = getDefaultParameters()
    % getDefaultParameters returns default tracking parameters.
    param.motionModel = 'ConstantAcceleration';
    param.initialLocation = 'Same as first detection';
    param.initialEstimateError = 1E5 * ones(1, 3);
    param.motionNoise = [25, 10, 1];
    param.measurementNoise = 25;
    param.segmentationThreshold = 0.05;
    param.entranceThreshold = 50;
    param.exitThreshold = 1230;
    param.outputFileName = 'trackedVideo_default.mp4';
    param.useKalmanFilter = true;
end
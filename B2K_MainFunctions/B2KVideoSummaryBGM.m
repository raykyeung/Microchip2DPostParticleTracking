function [L_start, R_end, tFrames, nFrames, ...
    flag_pChip_present, flag_pChip_stuck, flag_no_pChip, flag_video_corrupt, ...
    myFigures] = B2KVideoSummaryBGM(vidObj, myFigures, flagDebug, flagRunWithDisplay, addedFrames)
%% B2KVideoSummaryBGM.m - [Function] Batch single p-Chip tracking
%
% Author: Raymond Yeung
% Release date: 2025
% E-mail: ryeun003@ucr.edu
% B2K Group, Dept. of Bioengineering, Univ. of California, Riverside
% Victor G. J. Rodgers Dept. of Bioengineering, Univ. of California, Riverside
% William H. Grover, Dept. of Bioengineering, Univ. of California, Riverside
% Philip L. Brisk Dept. of Computer Science, Univ. of California, Riverside

%%
% Summarizes video frames using a foreground detector and blob analysis.
%
% This function processes a video to determine the frame range during which an object 
% (p‑Chip) is present. It first processes the video sequentially to obtain a detection 
% vector (ensuring proper training of the foreground detector), then uses a coarse-to-fine 
% search to quickly and robustly determine the first and last frames with detection.
%
% Only two flag values are used: 1 indicates true, and NaN indicates false or not applicable.
%
% OUTPUTS:
%   L_start            - Start frame index for analysis (with added margin).
%   R_end              - End frame index for analysis (with added margin).
%   tFrames            - Total number of frames in the summary.
%   nFrames            - Total number of frames in the video.
%   flag_pChip_present - 1 if p‑Chip is detected from the beginning then leaves; otherwise NaN.
%   flag_pChip_stuck   - 1 if p‑Chip enters later and persists until the end (or is present throughout); otherwise NaN.
%   flag_no_pChip      - 1 if no p‑Chip is detected in the video; otherwise NaN.
%   flag_video_corrupt - 1 if the video is corrupt (no frames); otherwise NaN.
%   myFigures          - Updated figure structure.

%% Initialize Flags and Check Video Validity
% (All flags are either 1 or NaN.)
flag_no_pChip      = NaN;
flag_pChip_present = NaN;
flag_pChip_stuck   = NaN;
flag_video_corrupt = NaN;

nFrames = vidObj.NumFrames;
if nFrames == 0
    flag_video_corrupt = 1;
    myFigures = B2KMyFigures(myFigures, flagRunWithDisplay);
    L_start = NaN;
    R_end   = NaN;
    tFrames = NaN;
    return;
else
    flag_video_corrupt = NaN;
end

%% Parameter Definitions
segmentationThreshold = 0.05;  % Fixed initial variance for the foreground detector.
coarseStep = 10;               % Coarse sampling step size.

%% Precompute Detection Results Sequentially
% Create a single detector and blob analyzer to process the video sequentially.
fgDetector = vision.ForegroundDetector('NumTrainingFrames', 10, ...
    'InitialVariance', segmentationThreshold);
blobAnalyzer = vision.BlobAnalysis('MinimumBlobArea', 50, 'MaximumBlobArea', 500, ...
    'PerimeterOutputPort', true, 'OrientationOutputPort', true, 'ExcludeBorderBlobs', true);

% Preallocate a logical vector to store detection (true if object is detected).
detectionFlag = false(nFrames, 1);

for i = 1:nFrames
    frame = read(vidObj, i);
    grayFrame = im2gray(im2single(frame));
    mask = step(fgDetector, grayFrame);
    mask = imopen(mask, strel('disk', 3));
    [area, ~, ~, ~, ~] = step(blobAnalyzer, mask);
    detectionFlag(i) = ~isempty(area);
    
    if flagDebug
        fprintf('Sequential processing: Frame %d, detection = %d\n', i, detectionFlag(i));
    end
end

%% Coarse Sampling of the Detection Vector
coarseIndices = (1:coarseStep:nFrames).'; % Frames to sample.
coarseDetection = detectionFlag(coarseIndices);

%% Identify Coarse Boundaries and Refine Them
if ~any(coarseDetection)
    % No object detected anywhere.
    flag_no_pChip = 1;
    L_start = NaN;
    R_end   = NaN;
    tFrames = NaN;
else
    flag_no_pChip = NaN;
    
    % Find the first and last indices (in the coarse sample) with detection.
    firstCoarseIdx = coarseIndices(find(coarseDetection, 1, 'first'));
    lastCoarseIdx  = coarseIndices(find(coarseDetection, 1, 'last'));
    
    % Refine the left boundary:
    % Search linearly from max(1, firstCoarseIdx - coarseStep + 1) to firstCoarseIdx.
    searchStart = max(1, firstCoarseIdx - coarseStep + 1);
    L_start_candidate = firstCoarseIdx;
    for i = searchStart:firstCoarseIdx
        if detectionFlag(i)
            L_start_candidate = i;
            if flagDebug
                fprintf('Fine refinement (left): Detection first found at frame %d\n', i);
            end
            break;
        end
    end
    
    % Refine the right boundary:
    % Search linearly from lastCoarseIdx to min(nFrames, lastCoarseIdx + coarseStep - 1).
    searchEnd = min(nFrames, lastCoarseIdx + coarseStep - 1);
    R_end_candidate = lastCoarseIdx;
    for i = lastCoarseIdx:searchEnd
        if detectionFlag(i)
            R_end_candidate = i;
            if flagDebug
                fprintf('Fine refinement (right): Detection confirmed at frame %d\n', i);
            end
        else
            if flagDebug
                fprintf('Fine refinement (right): Lost detection at frame %d\n', i);
            end
            break;
        end
    end
    
    % Apply the added margin.
    L_start = max(1, L_start_candidate - addedFrames);
    R_end   = min(nFrames, R_end_candidate + addedFrames);
    tFrames = R_end - L_start + 1;
    
    % Set flags based on the refined boundaries.
    if L_start_candidate == 1 && R_end_candidate < nFrames
        flag_pChip_present = 1;
        flag_pChip_stuck = NaN;
    elseif L_start_candidate > 1 && R_end_candidate == nFrames
        flag_pChip_stuck = 1;
        flag_pChip_present = NaN;
    elseif L_start_candidate == 1 && R_end_candidate == nFrames
        flag_pChip_stuck = 1;
        flag_pChip_present = NaN;
    else
        flag_pChip_present = 1;
        flag_pChip_stuck = NaN;
    end
end

%% Optional Display of Summary Figures
if flagRunWithDisplay
    myFigures = B2KMyFigures(myFigures, flagRunWithDisplay);
    figure(myFigures.fig.(myFigures.handle(myFigures.figNum)));
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot the detection vector with the refined boundaries.
    nexttile;
    plot(1:nFrames, double(detectionFlag), 'LineWidth', 1.5);
    hold on;
    yline(0.5, '--r', 'Refined Boundaries');
    xline(L_start, '--m', sprintf('Start: %d', L_start));
    xline(R_end, '--g', sprintf('End: %d', R_end));
    hold off;
    xlabel('Frame Number');
    ylabel('Detection Flag');
    title('Detection Across Frames');
    ylim([-0.1, 1.1]);
    
    if isnan(flag_no_pChip)
        % Display start frame.
        nexttile;
        frameStart = read(vidObj, L_start);
        imshow(frameStart);
        title(sprintf('Start Frame: %d', L_start));
        
        % Display a middle frame.
        nexttile;
        midFrameIdx = floor((L_start + R_end) / 2);
        frameMid = read(vidObj, midFrameIdx);
        imshow(frameMid);
        title(sprintf('Middle Frame: %d', midFrameIdx));
        
        % Display end frame.
        nexttile;
        frameEnd = read(vidObj, R_end);
        imshow(frameEnd);
        title(sprintf('End Frame: %d', R_end));
    else
        % Display start frame.
        nexttile;
        frameStart = read(vidObj, 1);
        imshow(frameStart);
        title(sprintf('Start Frame: %d', 1));
        
        % Display a middle frame.
        nexttile;
        midFrameIdx = floor((1 + nFrames) / 2);
        frameMid = read(vidObj, midFrameIdx);
        imshow(frameMid);
        title(sprintf('Middle Frame: %d', midFrameIdx));
        
        % Display end frame.
        nexttile;
        frameEnd = read(vidObj, nFrames);
        imshow(frameEnd);
        title(sprintf('End Frame: %d', nFrames));
    end
end

%% Clear figure helper if needed.
clear B2KMyFigures

end
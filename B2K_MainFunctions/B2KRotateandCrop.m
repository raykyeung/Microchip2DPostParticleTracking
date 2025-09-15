function [ CI, T ] = B2KRotateandCrop( I, ang )
%ROTATEANDCROP Rotate an image 'I' by 'ang' degrees, and crop its biggest
% inner rectangle.

[h,w,~] = size(I);
ang = deg2rad(ang);

% Affine rotation
R = [cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1];
T = affine2d(R);
B = imwarp(I,T);

% Largest rectangle
% solution from https://stackoverflow.com/a/16778797

wb = w >= h;
sl = w*wb + h*~wb;
ss = h*wb + w*~wb;

cosa = abs(cos(ang));
sina = abs(sin(ang));

if ss <= 2*sina*cosa*sl
    x = .5*min([w h]);
    wh = wb*[x/sina x/cosa] + ~wb*[x/cosa x/sina];
else
    cos2a = (cosa^2) - (sina^2);
    wh = [(w*cosa - h*sina)/cos2a (h*cosa - w*sina)/cos2a]; 
end

hw = flip(wh);

% Top-left corner
% tl = round(max(size(B)/2 - hw/2,1));
% tl = ceil(max(size(B)/2 - hw/2,1));
tl = ceil(size(B)/2 - hw/2) + 1;

% tly1 = size(B/2);
% tly2 = hw/2;
% 
% tly = round(max(tly1(1) - tly2(1)));
% tlx = round(max(tly1(2) - tly2(2)));
% tl = [tly tlx];

% Bottom-right corner
% br = tl + round(hw);
br = tl + floor(hw) - 2;

% Cropped image
CI = B(tl(1):br(1),tl(2):br(2),:);
% CI = B(tl(1)+2:br(1),tl(2):br(2),:);
% CI = B(tl(1)+2:br(1),tl(2)-2:br(2)+2,:);
% CI = B(tl(1)+2:br(1),tl(2)-4:br(2)+4,:);
% CI = B(tl(1)+2:br(1),tl(2)-8:br(2)+8,:);
% CI = B(tl(1)+2:br(1),tl(2)-50:br(2)+50,:);
% CI = B(tl(1)+2:br(1),tl(2)-30:br(2)+30,:);
% CI = B(tl(1)+2:br(1),tl(2)-14:br(2)+14,:);
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
tl = ceil(size(B)/2 - hw/2) + 1;

% Bottom-right corner
br = tl + floor(hw) - 2;

% Cropped image
CI = B(tl(1):br(1),tl(2):br(2),:);

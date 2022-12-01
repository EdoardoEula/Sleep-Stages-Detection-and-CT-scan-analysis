clc
clear
close all

%% All times extraction
t = [{'0'}, {'10'}, {'20'}, {'30'}, {'40'}, {'50'}, {'60'}, {'70'}, {'80'}, {'90'}];
vol = cell(length(t), 1);
info = cell(length(t), 1);
for ti = 1:length(t)
    [vol{ti}, info{ti}] = readDCMfolder(t{ti});
end
%% Example image
im = vol{1}(:,:,1);
%%
I_max = max(max(im));
im_bin = imbinarize(im, I_max/2);

figure
C = imcontour(im, 1, 'b');

figure
subplot(1,3,1), imshow(im_bin,[]), title('Binarized')
subplot(1,3,2), imcontour(im,1,'b')
subplot(1,3,3), imshow(im_bin,[]), hold on, imcontour(im,1,'m'), title('Contours')

%%
figure
imshow(im, [])
h = drawellipse('Center',[267 250],'SemiAxes',[78 72], ...
    'RotationAngle',287,'StripeColor','m');
%%
mask = createMask(h);
figure
imshow(mask)

%%
I_max = max(max(im));
im_bin = imbinarize(im, I_max/2);
im_bin_mask = im_bin .* mask;

figure
imshow(im_bin_mask, [])
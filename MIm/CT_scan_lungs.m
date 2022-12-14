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
% I_max = max(max(im));
% im_bin = imbinarize(im, I_max/2);
% 
% figure
% C = imcontour(im, 1, 'b');
% 
% figure
% subplot(1,3,1), imshow(im_bin,[]), title('Binarized')
% subplot(1,3,2), imcontour(im,1,'b')
% subplot(1,3,3), imshow(im_bin,[]), hold on, imcontour(im,1,'m'), title('Contours')

%%
figure
imshow(im, [])
[X,Y] = getpts();
% h = drawellipse('Center',[260 267],'SemiAxes',[140 220], ...
%     'RotationAngle',90,'StripeColor','m'); %sistemare misure ellisse

%%
mask = poly2mask(X, Y, 512, 512);

%%
I_max = max(max(im));
%im_bin = imbinarize(im, I_max/2);
%im_bin_mask = im_bin .* mask;
roi = im .* mask;

figure
subplot(121), imshow(im, []);
subplot(122), imshow(roi, []);

%%
%Normalisation
max_im = max(im(:));
min_im = min(im(:));
im = (im-min_im)./(max_im-min_im);
figure
imshow(im, []);

figure
subplot(311), imhist(im, 32), title('N = 32');
subplot(312), imhist(im, 64), title('N = 64');
subplot(313), imhist(im, 256), title('N = 256');

pause(3)
close

%%
% gamma_vect = 0.7:0.1:1.5;
LOW_IN = min(im(:));
HIGH_IN = max(im(:));

LOW_OUT = 0;
HIGH_OUT = 1;
% 
% roi_adj = zeros(size(roi, 1), size(roi, 2), 1, length(gamma_vect));
% 
% for ind = 1:length(gamma_vect)
%     roi_adj(:, :, 1, ind) = imadjust(roi, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma_vect(ind));
% end
% 
% figure, montage(roi_adj)

%% Adjusted image
gamma = 0.9;
roi_adj = imadjust(im, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma);
figure, 
subplot(121), imshow(roi_adj)
subplot(122), imhist(roi_adj)

%% T-transformation
roi_adj_t = roi_adj;
roi_adj_t(roi_adj > 0.65) = 0;

figure, 
subplot(121), imshow(roi_adj_t)
subplot(122), imhist(roi_adj_t)

%%
figure, imshow(roi_adj_t)
[X,Y] = getpts()
bw = bwselect(roi_adj_t, X, Y)

%%
figure, imshow(bw)

%%
segmented_image(512,512,100) = 0;
bw_t(512,512,59) = 0;
for i = 1:59
    im = vol{1}(:,:,i);
    max_im = max(im(:));
    min_im = min(im(:));
    im = (im-min_im)./(max_im-min_im);
    XY = im;
    gamma = 0.8;
    im=XY;
    LOW_IN = min(im(:));
    HIGH_IN = max(im(:));
    LOW_OUT = 0;
    HIGH_OUT = 1;
    roi_adj = imadjust(im, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma);
    roi_adj_t = roi_adj;
    roi_adj_t(roi_adj > 0.65) = 0;
    bw = bwselect(roi_adj_t, X, Y);
    bw = imfill(bw,"holes");
    figure, imshow(im.*single(bw))
    %roi_adj_tN(:,:,i) = imcomplement(roi_adj_t);
    im = histeq(im);
    segmented_image(:,:,i) = im.*single(bw);
    bw_t(:,:,i) = im.*single(bw);
end

%%
volumeViewer(segmented_image)

%%

volLungsPixels = regionprops3(bw_t,"volume");

spacingx = 0.76;
spacingy = 0.76;
spacingz = 1.26*1e-6;
unitvol = spacingx*spacingy*spacingz;

volLungs = volLungsPixels.Volume(1)*unitvol;

%% Cross-sectional area
V = vol{1}(:,:,:);

%Normalization
for i=1:136
    max_im = max(V(:,:,i));
    min_im = min(V(:,:,i));
    V(:,:,i) = (V(:,:,i)-min_im)./(max_im-min_im);
end

%We choose the middle slice
XY = V(:,:,30);
XZ = squeeze(V(256,:,:));

imshow(XZ, [])

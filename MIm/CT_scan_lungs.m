clc
clear
close all

%% All times extraction
t = [{'0'}, {'10'}, {'20'}, {'30'}, {'40'}, {'50'}, {'60'}, {'70'}, {'80'}, {'90'}];
vol = cell(length(t), 1);
info = cell(length(t), 1);
for ti = 1:length(t)
    [vol{ti}, ~] = readDCMfolder(t{ti});
end

%% Dicom info
t = '0';
[~, info] = readDCMfolder(t);
%% t = 0
t = 1;
im = vol{t}(:,:,1);
%Normalisation
max_im = max(im(:));
min_im = min(im(:));
im = (im-min_im)./(max_im-min_im);
figure
[im, RECT] = imcrop(im);

%%
gamma_vect = 0.6:0.1:1.5;
LOW_IN = min(im(:));
HIGH_IN = max(im(:));
LOW_OUT = 0;
HIGH_OUT = 1;
im_adj= zeros(size(im,1),size(im,2),length(gamma_vect));
for i = 1:length(gamma_vect)
    im_adj(:,:,i) = imadjust(im, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma_vect(i));
end

figure, montage(im_adj)
%%
gamma = 0.6;

im_adj = imadjust(im, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma);
figure, imhist(im_adj)
figure, imshow(im_adj, [], 'InitialMagnification', 'fit')
[X,Y] = getpts();

%% Segmentation over time = 0
segmented_image = zeros(size(im,1),size(im,2),59);
bw_t = zeros(size(im,1),size(im,2),59);
for i = 1:59
    im = vol{t}(:,:,i);
    im = imcrop(im, RECT);
    max_im = max(im(:));
    min_im = min(im(:));
    im = (im-min_im)./(max_im-min_im);
    LOW_IN = min(im(:));
    HIGH_IN = max(im(:));
    LOW_OUT = 0;
    HIGH_OUT = 1;
    im_adj = imadjust(im, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma);

    im_adj(im_adj > 0.5) = 0;
    bw = bwselect(im_adj, X, Y);
    bw = imfill(bw,"holes");
    %figure, imshow(imcomplement(im.*single(bw)), [], 'InitialMagnification', 'fit')
    segmented_image(:,:,i) = im.*single(bw);
    bw_t(:,:,i) = single(bw);
end

%%
for i = 1:55
    slice = segmented_image(:,:,i);
    AreaLungsPixels = regionprops(logical(slice), 'Area');
    AreaLungsPixels = struct2table(AreaLungsPixels);
    spacingx = 0.976;
    spacingy = 0.976;
    spacingz = 2.5;
    unitar = spacingx*spacingy;
    unitvol = spacingx*spacingy*spacingz;
    AreaSx(i) = (AreaLungsPixels.Area(2)*unitar)*1e-2;
    AreaDx(i) = (AreaLungsPixels.Area(1)*unitar)*1e-2;
end

figure, plot(1:55, AreaSx, 1:55, AreaDx)
legend('Sx', 'Dx')
title('Area')

%%
volumeViewer(segmented_image)

%%
volLungsPixels = regionprops3(logical(bw_t),"volume");

spacingx = 0.976;
spacingy = 0.976;
spacingz = 2.5*1e-6;
unitvol = spacingx*spacingy*spacingz;

volLungs = volLungsPixels.Volume(1)*unitvol + volLungsPixels.Volume(2)*unitvol;


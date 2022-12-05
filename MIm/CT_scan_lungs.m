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
im = vol{1}(:,:,50);
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
h = drawellipse('Center',[260 267],'SemiAxes',[140 220], ...
    'RotationAngle',90,'StripeColor','m'); %sistemare misure ellisse
%%
mask = createMask(h);

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
max_im = max(roi(:));
min_im = min(roi(:));
roi = (roi-min_im)./max_im;
figure
imshow(roi, []);

figure
subplot(311), imhist(roi, 32), title('N = 32');
subplot(312), imhist(roi, 64), title('N = 64');
subplot(313), imhist(roi, 256), title('N = 256');

pause(3)
close

%%
% gamma_vect = 0.7:0.1:1.5;
% LOW_IN = min(roi(:));
% HIGH_IN = max(roi(:));
% 
% LOW_OUT = 0;
% HIGH_OUT = 1;
% 
% roi_adj = zeros(size(roi, 1), size(roi, 2), 1, length(gamma_vect));
% 
% for ind = 1:length(gamma_vect)
%     roi_adj(:, :, 1, ind) = imadjust(roi, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma_vect(ind));
% end
% 
% figure, montage(roi_adj)

%% Adjusted image
gamma = 0.3;
roi_adj = imadjust(roi, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma);
figure, 
subplot(121), imshow(roi_adj)
subplot(122), imhist(roi_adj)

%% T-transformation
roi_adj_t = roi_adj;
roi_adj_t(roi_adj > 0.65) = 0;

roi_adj_tN = imcomplement(roi_adj_t);

figure, 
subplot(121), imshow(roi_adj_tN)
subplot(122), imhist(roi_adj_t)

%%
roi_adj_tN(:,:,136) = 0;
for i = 1:136
    im = vol{1}(:,:,i);
    h = drawellipse('Center',[260 267],'SemiAxes',[140 220], ...
        'RotationAngle',90,'StripeColor','m');
    mask = createMask(h);
    I_max = max(max(im));
    roi = im .* mask;
    max_im = max(roi(:));
    min_im = min(roi(:));
    roi = (roi-min_im)./max_im;
    gamma = 0.3;
    roi_adj = imadjust(roi, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma);
    roi_adj_t = roi_adj;
    roi_adj_t(roi_adj > 0.65) = 0;
    roi_adj_tN(:,:,i) = imcomplement(roi_adj_t);
end

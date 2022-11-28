clc
clear
close all

%% All times extraction
t = [{'0'}, {'10'}, {'20'}, {'30'}, {'40'}, {'50'}, {'60'}, {'70'}, {'80'}, {'90'}];
D = cell(length(t), 1);
vol = cell(length(t), 1);
info = cell(length(t), 1);
for ti = 1:length(t)
    [vol{ti}, info{ti}] = readDCMfolder(t{ti});
end

%% Example image
im = vol{1}(:,:,1);

%%
figure
imshow(im, [])
[X,Y] = getpts;

%%



=======
%% T0
t = '0';
[D, vol,info] = readDCMfolder(t);
>>>>>>> parent of 12c78e6 (Files)

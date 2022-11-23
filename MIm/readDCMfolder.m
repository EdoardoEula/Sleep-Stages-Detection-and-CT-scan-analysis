function [D, vol,info] = readDCMfolder(t)
%%
folder = ['Patient0/T_', t, '/CT/'];
D = dir([folder, '*.dcm']);

im = dicomread([folder, D(1).name]);
vol = zeros(size(im,1), size(im,2), size(D,1));

info = dicominfo([folder, D(1).name]);

f = waitbar(0, sprintf('Loading:  %u / %u', 0, size(D,1)));
for ind=1:size(D, 1)
    waitbar(ind/size(D,1), f, sprintf('Loading:  %u / %u', ind, size(D,1)));
    vol(:,:,ind) = dicomread( [folder, D(ind).name] );
end
close(f);

end
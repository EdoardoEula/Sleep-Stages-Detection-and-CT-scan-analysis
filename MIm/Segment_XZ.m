function [segmented_image_XZ] = Segment_XZ(vol,T, gamma, X, Z)
for t = 1:length(T)
    if t~=6
        im_vect_XZ = vol{t}(180:315,:,1:52);
    else 
        im_vect_XZ = vol{t}(180:315,:, 85:136);
    end 

    for i = 1:136
        im_XZ = squeeze(im_vect_XZ(i,:,:));
        max_im = max(im_XZ(:));
        min_im = min(im_XZ(:));
        im_XZ = (im_XZ-min_im)./(max_im-min_im);
        LOW_IN = min(im_XZ(:));
        HIGH_IN = max(im_XZ(:));
        LOW_OUT = 0;
        HIGH_OUT = 1;
        im_adj_XZ = imadjust(im_XZ, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma);
    
        im_adj_XZ(im_adj_XZ > 0.55) = 0;
        bw_XZ = bwselect(im_adj_XZ, X, Z);
        bw_XZ = imfill(bw_XZ,"holes"); %correction of holes
        segmented_image_XZ{t}(i,:,:) = im_XZ.*single(bw_XZ);
    end 
end
end


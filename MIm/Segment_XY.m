function [segmented_image] = Segment_XY(vol,T, gamma, RECT, X, Y, noise)
 for t=1:length(T)
    if t~=6
        im_vect = vol{t}(:,:,1:52);
    else 
        im_vect = vol{t}(:,:, 85:136);
    end 
    for i = 1:52
        im = im_vect(:,:,i);
        im = imcrop(im, RECT);
        max_im = max(im(:));
        min_im = min(im(:));
        im = (im-min_im)./(max_im-min_im);
        switch noise
            case 1
                im = imnoise(im, 'gaussian', 0, 1e-3);
            case 2
                im = imnoise(im, 'gaussian', 0, 0.1);
            case 3
                im = imnoise(im, 'gaussian', 0.5, 1e-3);
            case 4
                im = imnoise(im, 'salt & pepper', 0.05);
            case 5
                im = imnoise(im, 'salt & pepper', 0.2);
            case 6
                im = imnoise(im, 'salt & pepper', 0.35);
            otherwise
                im = im;
        end

        LOW_IN = min(im(:));
        HIGH_IN = max(im(:));
        LOW_OUT = 0;
        HIGH_OUT = 1;
        im_adj = imadjust(im, [LOW_IN HIGH_IN], [LOW_OUT HIGH_OUT], gamma);
    
        im_adj(im_adj > 0.5) = 0;
        bw = bwselect(im_adj, X, Y);
        bw = imfill(bw,"holes"); %correction of holes
        segmented_image{t}(:,:,i) = im.*single(bw);
    end
end 
end


function [Area_tot] = Area_XY(segmented_image,T)

for t=1:length(T)
    for i = 1:52
        slice = segmented_image{t}(:,:,i);
        AreaLungsPixels = regionprops(logical(slice), 'Area');
        AreaLungsPixels = struct2table(AreaLungsPixels);
        spacingx = 0.976;
        spacingy = 0.976;
        unitar = spacingx*spacingy;
        Area_tot{t}(i)=(sum(AreaLungsPixels.Area)*unitar)*1e-2;
    end
end
end


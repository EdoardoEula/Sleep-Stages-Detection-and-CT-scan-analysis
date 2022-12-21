function [volLungs] = Volume(segmented_image, T)
for t = 1:length(T)

    volLungsPixels = regionprops3(logical(segmented_image{t}),"volume");

    spacingx = 0.976;
    spacingy = 0.976;
    spacingz = 2.5*1e-6;
    unitvol = spacingx*spacingy*spacingz;

    volLungs(t) = volLungsPixels.Volume(1)*unitvol + volLungsPixels.Volume(2)*unitvol;
end
end


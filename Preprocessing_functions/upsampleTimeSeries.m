dicoms = double(readDICOMS);

desiredFrames = 40;
resampled = zeros(size(dicoms,1),size(dicoms,2),desiredFrames);
for i=1:size(dicoms,1)
    for j=1:size(dicoms,2)
        pixelLine = squeeze(dicoms(i,j,:));
        x = 1:length(pixelLine);
        xi = linspace(1,length(pixelLine),desiredFrames);
        yi = interp1(x,pixelLine,xi,'spline');
        resampled(i,j,:) = yi;
    end 
end 

        
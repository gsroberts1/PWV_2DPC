% Load data, enter wrapped frames, set Venc
loadRadial2DPC;
frames = 20:29;
Venc = 1500;

% Draw ROI for all frames
figure; imshow(MAG,[]);
circle = drawcircle;
radius = circle.Radius; 
center = round(circle.Center);
[X,Y] = ndgrid(1:size(MAG,1),1:size(MAG,2));
X = X-center(2); %shift coordinate grid
Y = Y-center(1);
roiMask = sqrt(X.^2+Y.^2)<=radius;

% Unwrap each frame
for i=1:length(frames)
    frame = frames(i);
    plane = vz(:,:,frame);
    
    newPlane = plane.*roiMask;
    newPlane = (newPlane/Venc)*pi;
    UW = Unwrap_TIE_DCT_Iter(newPlane);
    UW = (UW/pi)*Venc;
    %UW = UW-3000;
    %UW(UW==-3000) = 0;
    %UW = newROI.*UW;
    %UW = padarray(UW,[226 226]);
    figure; imshowpair(newPlane,UW,'montage');

    plane = plane.*(~roiMask);
    plane = plane + UW;

    % Overwrite existing dat file
    filename = ['ph_' num2str(frame-1,'%03.f') '_vd_3.dat'];
    fid = fopen(filename,'w');
    fwrite(fid,plane,'int16');
    fclose(fid);
end 
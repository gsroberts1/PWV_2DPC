%% Data Loading
directory = pwd;
fid = fopen([directory '\pcvipr_header.txt'], 'r');
dataArray = textscan(fid, '%s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fid); clear ans;

dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);
resx = pcviprHeader.matrixx;  
resy = pcviprHeader.matrixy;  
nframes = pcviprHeader.frames;           

vz  = zeros(resx,resy,nframes);
for j = 1:nframes   
    vz(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_vd_3.dat']),[resx resy]);
end
vz = squeeze(vz);

%% Unwrapping
frames = 16:20;
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
    UW = (UW/pi)*-Venc;
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
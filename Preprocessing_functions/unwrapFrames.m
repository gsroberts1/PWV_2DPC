clean;
load2DPC;
frame = 18;
plane = v(:,:,3,frame);

figure; imshow(MAG,[]);
circle = drawcircle;
radius = circle.Radius; 
center = round(circle.Center);

[X,Y] = ndgrid(1:size(plane,1),1:size(plane,2));
X = X-center(2); %shift coordinate grid
Y = Y-center(1);
roiMask = sqrt(X.^2+Y.^2)<=radius;
newPlane = plane.*roiMask;
%newPlane = newPlane(center(2)-30:center(2)+29,center(1)-30:center(1)+29);
newPlane = (newPlane/1500)*pi;
%newROI = roiMask(center(2)-30:center(2)+29,center(1)-30:center(1)+29);
UW = Unwrap_TIE_DCT_Iter(newPlane);
UW = (UW/pi)*1500;
%UW = UW-3000;
%UW(UW==-3000) = 0;
%UW = newROI.*UW;
%UW = padarray(UW,[226 226]);
figure; imshowpair(newPlane,UW,'montage');

plane = plane.*(~roiMask);
plane = plane + UW;

figure; imshow(plane,[]);
filename = ['ph_' num2str(frame-1,'%03.f') '_vd_3.dat'];
fid = fopen(filename,'w');
fwrite(fid,plane,'int16');
fclose(fid);
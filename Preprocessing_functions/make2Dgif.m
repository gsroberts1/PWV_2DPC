%load2DPC;
%vz = squeeze(v(:,:,3,:));
% figure; imshow(fiestaAAo(:,:,1),[]);
% 
% circle = drawcircle;
% radius = 25; %get radius of circle
% center = round(circle.Center); %get center coordinates
% [X,Y] = ndgrid(1:size(fiestaAAo,1),1:size(fiestaAAo,2));
% X = X-center(2); %shift coordinate grid
% Y = Y-center(1);
% roiMask = sqrt(X.^2+Y.^2)<=radius; %anything outside radius is ignored
% 
% circle2 = drawcircle;
% radius2 = 25; %get radius of circle
% center2 = round(circle2.Center); %get center coordinates
% [X,Y] = ndgrid(1:size(fiestaAAo,1),1:size(fiestaAAo,2));
% X = X-center2(2); %shift coordinate grid
% Y = Y-center2(1);
% roiMask2 = sqrt(X.^2+Y.^2)<=radius2; %anything outside radius is ignored
% 
% rois = roiMask + roiMask2;
%close all;

test = pwvAAo + 1000;
test = test./2000;
test(test>1) = 1;
test(test<0) = 0;
test = uint8(test.*255);

for i=1:size(pwvAAo,3)
    s1 = regiongrowing(fiestaAAo(:,:,i),212,212,240);
    s1 = imclose(s1,strel('octagon',12));
    s1 = imopen(s1,strel('disk',5));
    s1 = imdilate(s1,strel('disk',1));
    s2 = regiongrowing(fiestaAAo(:,:,i),314,268,165);
    s2 = imclose(s2,strel('disk',7));
    s2 = imopen(s2,strel('disk',7));
    s2 = imdilate(s2,strel('disk',2));
    if i>35
        figure; imshow(s2);
        circle = drawcircle;
        radius = circle.Radius; %get radius of circle
        center = round(circle.Center); %get center coordinates
        [X,Y] = ndgrid(1:size(fiestaAAo,1),1:size(fiestaAAo,2));
        X = X-center(2); %shift coordinate grid
        Y = Y-center(1);
        roiMask = sqrt(X.^2+Y.^2)<=radius; %anything outside radius is ignored

        
        s2 = roiMask + s2;
    end 
    rois = s1 + s2; 
    

    
    y = colormap(jet);
    imwrite(test(:,:,i),y,'rgb.jpeg','Mode','lossless'); % Re-change it to colored one
    colored = imread('rgb.jpeg');
    imshow(fiestaAAo(:,:,i),[]); colormap('gray'); hold on; imagesc(colored,'AlphaData',rois);
    F(i) = getframe();
    close all;
end 




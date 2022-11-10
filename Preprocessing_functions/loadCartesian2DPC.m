dirinfo = [dir('*.dcm'); dir('*.dicom'); dir('*.sdcopen')];
if isempty(dirinfo)
    disp('NO DICOMS FOUND IN DIRECTORY');
end 

header = dicominfo(dirinfo(1).name);
for file = 1:length(dirinfo)
    DICOMS(:,:,file) = dicomread(dirinfo(file).name);
end 
mag = DICOMS(:,:,41:end);
vz = DICOMS(:,:,1:40);

originShift = [header.ImagePositionPatient;1]; % origin is top left corner of image
xres = header.PixelSpacing(1);
yres = header.PixelSpacing(2);
zres = header.SliceThickness;

% sometimes get extremely small values that should be 0, so round
xVector = round(header.ImageOrientationPatient(1:3),8); % what direction rows run w/r/to x
yVector = round(header.ImageOrientationPatient(4:6),8); % what direction the cols run w/r/to y
zVector = [cross(xVector,yVector);0];

xVector = [xVector;0];
yVector = [yVector;0];
rotation = [xres*xVector yres*yVector zres*zVector originShift]; % turn vectors into matrices

spatialRes = header.PixelSpacing(1);
temporalRes = header.NominalInterval/header.CardiacNumberOfImages;
disp(['Spatial Resolution = ' num2str(spatialRes) ' mm']);
disp(['Temporal Resolution = ' num2str(temporalRes) ' ms']);
clear dirinfo file header 
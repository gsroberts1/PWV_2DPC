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

spatialRes = header.PixelSpacing(1);
temporalRes = header.NominalInterval/header.CardiacNumberOfImages;
disp(['Spatial Resolution = ' num2str(spatialRes) ' mm']);
disp(['Temporal Resolution = ' num2str(temporalRes) ' ms']);
clear dirinfo file header 
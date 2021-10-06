directory = pwd;
fid = fopen([directory '\pcvipr_header.txt'], 'r');
dataArray = textscan(fid, '%s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fid); clear ans;

dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);
resx = pcviprHeader.matrixx;  
resy = pcviprHeader.matrixy;  
nframes = pcviprHeader.frames;           

MAG = load_dat(fullfile(directory,'MAG.dat'),[resx resy]);
CD = load_dat(fullfile(directory,'CD.dat'),[resx resy]);
VMEAN = load_dat(fullfile(directory,'comp_vd_3.dat'),[resx resy]);

vz  = zeros(resx,resy,nframes);
mag = zeros(resx,resy,nframes);
cd  = zeros(resx,resy,nframes);
for j = 1:nframes   
    vz(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_vd_3.dat']),[resx resy]);
    mag(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_mag.dat']),[resx resy]);
    cd(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_cd.dat']),[resx resy]);
end
MAG = flipud(MAG);
CD = flipud(CD);
VMEAN = flipud(VMEAN);
mag = flipud(mag);
cd = flipud(cd);
vz = flipud(squeeze(vz));

spatialRes = nonzeros(abs(([pcviprHeader.ix pcviprHeader.iy pcviprHeader.iz])));
temporalRes = pcviprHeader.timeres;
disp(['Spatial Resolution = ' num2str(spatialRes) ' mm']);
disp(['Temporal Resolution = ' num2str(temporalRes) ' ms']);
clear dataArray directory fid j nframes pcviprHeader resx resy
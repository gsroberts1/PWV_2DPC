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

sx = pcviprHeader.sx;
sy = pcviprHeader.sy;
sz = pcviprHeader.sz;
originShift = [sx;sy;sz;1];
ix = pcviprHeader.ix;
iy = pcviprHeader.iy;
iz = pcviprHeader.iz;
jx = pcviprHeader.jx;
jy = pcviprHeader.jy;
jz = pcviprHeader.jz;
kx = pcviprHeader.kx;
ky = pcviprHeader.ky;
kz = pcviprHeader.kz;

xVector = round([ix;iy;iz;0],8); % what direction rows run w/r/to x
yVector = round([jx;jy;jz;0],8); % what direction the cols run w/r/to y
zVector = round([kx;ky;kz;0],8); % what direction the cols run w/r/to y
rotation = [xVector yVector zVector originShift];

spatialRes = nonzeros(abs(([pcviprHeader.ix pcviprHeader.iy pcviprHeader.iz])));
temporalRes = pcviprHeader.timeres;
disp(['Spatial Resolution = ' num2str(spatialRes) ' mm']);
disp(['Temporal Resolution = ' num2str(temporalRes) ' ms']);
clear dataArray directory fid j nframes pcviprHeader resx resy
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
VMEAN = zeros(resx,resy,3);
for n = 1:3
    VMEAN(:,:,n) = load_dat(fullfile(directory,['comp_vd_' num2str(n) '.dat']),[resx resy]);
end

v = zeros(resx,resy,3,nframes);
mag = zeros(resx,resy,nframes);
cd = zeros(resx,resy,nframes);
for j = 1:nframes   
    v(:,:,1,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_vd_1.dat']),[resx resy]);
    v(:,:,2,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_vd_2.dat']),[resx resy]);
    v(:,:,3,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_vd_3.dat']),[resx resy]);
    mag(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_mag.dat']),[resx resy]);
    cd(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_cd.dat']),[resx resy]);
end
vz = squeeze(v(:,:,3,:));

MAG = flipud(MAG);
CD = flipud(CD);
VMEAN = flipud(VMEAN);
v = flipud(v);
vz = flipud(vz);
mag = flipud(mag);
cd = flipud(cd);

spatialRes = nonzeros(abs(([pcviprHeader.ix pcviprHeader.iy pcviprHeader.iz])));
temporalRes = pcviprHeader.timeres;
disp(['Spatial Resolution = ' num2str(spatialRes) ' mm']);
disp(['Temporal Resolution = ' num2str(temporalRes) ' ms']);

function v = load_dat(name, res)
%% Load Dat
% Loads in dat files in current directory.
[fid,errmsg]= fopen(name,'r');
if fid < 0  % If name does not exist in directory
    disp(['Error Opening Data : ',errmsg]);
end

% Reads in as short, reshapes by image res.
v = reshape(fread(fid,'short=>single'),res);
fclose(fid);
end 
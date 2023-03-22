%% Check Reconstructed Image
directory = pwd;
fid = fopen([directory '\pcvipr_header.txt'], 'r');
dataArray = textscan(fid, '%s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fid); clear ans;

dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);
resx = pcviprHeader.matrixx;  
resy = pcviprHeader.matrixy;  
nframes = pcviprHeader.frames;        
spatialRes = nonzeros(abs(([pcviprHeader.ix pcviprHeader.iy pcviprHeader.iz])));
temporalRes = pcviprHeader.timeres;

MAG = load_dat(fullfile(directory,'MAG.dat'),[resx resy]);
CD = load_dat(fullfile(directory,'CD.dat'),[resx resy]);
VX = load_dat(fullfile(directory,'comp_vd_1.dat'),[resx resy]);
VY = load_dat(fullfile(directory,'comp_vd_2.dat'),[resx resy]);
VZ = load_dat(fullfile(directory,'comp_vd_3.dat'),[resx resy]);

vz  = zeros(resx,resy,nframes);
mag = zeros(resx,resy,nframes);
cd  = zeros(resx,resy,nframes);
for j = 1:nframes   
    vx(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_vd_1.dat']),[resx resy]);
    vy(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_vd_2.dat']),[resx resy]);
    vz(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_vd_3.dat']),[resx resy]);
    mag(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_mag.dat']),[resx resy]);
    cd(:,:,j) = load_dat(fullfile(directory, ['\ph_' num2str(j-1,'%03i') '_cd.dat']),[resx resy]);
end
MAG = flipud(MAG);
CD = flipud(CD);
VMEAN = flipud(VZ);
mag = flipud(mag);
cd = flipud(cd);
vz = flipud(squeeze(vz));

% sx = pcviprHeader.sx;
% sy = pcviprHeader.sy;
% sz = pcviprHeader.sz;
% originShift = [sx;sy;sz;1];
% ix = pcviprHeader.ix;
% iy = pcviprHeader.iy;
% iz = pcviprHeader.iz;
% jx = pcviprHeader.jx;
% jy = pcviprHeader.jy;
% jz = pcviprHeader.jz;
% kx = pcviprHeader.kx;
% ky = pcviprHeader.ky;
% kz = pcviprHeader.kz;
% 
% xVector = round([ix;iy;iz;0],8); % what direction rows run w/r/to x
% yVector = round([jx;jy;jz;0],8); % what direction the cols run w/r/to y
% zVector = round([kx;ky;kz;0],8); % what direction the cols run w/r/to y
% rotation = [xVector yVector zVector originShift];
figure; imshow3D(vz); 


%% Check Respiratory Binning
time = importdata('Time.txt');
resp = importdata('Resp.txt');
timeweight = importdata('TimeWeight.txt');
weight = importdata('Weight.txt');

last = round(length(time)/2);
time = time(1:last);
resp = resp(1:last);
weight = weight(1:last);

%window = 50;
%RESP = smoothdata(resp,'gaussian',window);
RESP = resp;
INSP = weight.*RESP;
INSP(INSP==0) = NaN;
EXP = (abs(weight-1)).*RESP;
EXP(EXP==0) = NaN;
pctDataUsed = sum(weight)/last*100;

figure; plot(time(1:last),RESP(1:last),'Color',[0.4940 0.1840 0.5560],'LineWidth',1.8); ...
    hold on; plot(time(1:last),EXP(1:last),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.8); ...
    hold on; plot(time(1:last),INSP(1:last),'Color',[0.4940 0.1840 0.5560],'LineWidth',1.8);
xlabel('Time (s)');
ylabel('Respiratory Amplitude (a.u.)');
title('Respiratory Waveform');
clear time* resp weight last EXP INSP RESP j fid dataArray


%% Command Window Output
disp(['Spatial Resolution = ' num2str(spatialRes) ' mm']);
disp(['Temporal Resolution = ' num2str(temporalRes) ' ms']);
disp(['Percent Data Used with Resp. Gating = ' num2str(round(pctDataUsed)) '%']);
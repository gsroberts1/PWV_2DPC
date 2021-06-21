function varargout = DrawCenterlines(varargin)
% DRAWCENTERLINES MATLAB code for DrawCenterlines.fig
%      DRAWCENTERLINES, by itself, creates a new DRAWCENTERLINES or raises the existing
%      singleton*.
%
%      H = DRAWCENTERLINES returns the handle to a new DRAWCENTERLINES or the handle to
%      the existing singleton*.
%
%      DRAWCENTERLINES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRAWCENTERLINES.M with the given input arguments.
%
%      DRAWCENTERLINES('Property','Value',...) creates a new DRAWCENTERLINES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DrawCenterlines_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DrawCenterlines_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DrawCenterlines

% Last Modified by GUIDE v2.5 28-May-2021 16:45:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DrawCenterlines_OpeningFcn, ...
                   'gui_OutputFcn',  @DrawCenterlines_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before DrawCenterlines is made visible.
function DrawCenterlines_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

% Update handles structure
set(handles.ImageSlider,'Enable','off'); 
set(handles.MinContrastUpdate,'Enable','off'); 
set(handles.MaxContrastUpdate,'Enable','off'); 
set(handles.WhiteBoxReminder,'Visible','off');
set(handles.ShowPlanesRadio,'Enable','off');
handles.CurrView = 'Axial';
handles.CurrCartesian = 1;
handles.CurrRadial = 1;

guidata(hObject, handles);
% UIWAIT makes DrawCenterlines wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DrawCenterlines_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in LoadCenterlinePush.
function LoadCenterlinePush_Callback(hObject, eventdata, handles)
set(handles.ImageSlider,'Enable','on'); 
set(handles.MinContrastUpdate,'Enable','on'); 
set(handles.MaxContrastUpdate,'Enable','on'); 
set(handles.WhiteBoxReminder,'Visible','off');

[clFile, clDir] = uigetfile({'*.mat','Useable Files (*.mat)';
   '*.mat',  'MATLAB files (*.mat)'; ...
   '*.*',  'All Files (*.*)'}, 'Select "anatCLdataset.mat" file');
load([clDir clFile]);
handles.axial = anatCLdataset.Axial;
handles.coronal = anatCLdataset.Coronal;
handles.sagittal = anatCLdataset.Sagittal;
handles.cartesian = anatCLdataset.Cartesian2DPC;
handles.Centerline = anatCLdataset.Centerline;
updateAnatImages(handles)

POINTS = handles.axial.POINTS;
axes(handles.CenterlineDisplay); hold on;
h = scatter3(POINTS(1,:),POINTS(2,:),POINTS(3,:), 32, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75], ...
    'DisplayName','Axial');
alpha = 0.3;
set(h,'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
legend('Location','southeast');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');

POINTS = handles.coronal.POINTS;
h = scatter3(POINTS(1,:),POINTS(2,:),POINTS(3,:), 32, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.75 .75 0], ...
    'DisplayName','Coronal');
alpha = 0.3;
set(h,'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);

POINTS = handles.sagittal.POINTS;
h = scatter3(POINTS(1,:),POINTS(2,:),POINTS(3,:), 32, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.75 0 .75], ...
    'DisplayName','Sagittal');
alpha = 0.3;
set(h,'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
legend('Location','southeast');

for t=1:length(handles.cartesian)
    POINTS = handles.cartesian(t).POINTS;
    scatter3(POINTS(1,:),POINTS(2,:),POINTS(3,:),'k*', ...
        'LineWidth',12, ...
        'DisplayName','Cart-2DPC');
    legend('Location','southeast','AutoUpdate','off');
end 

splineLine = handles.Centerline;
axes(handles.CenterlineDisplay);
plot3(splineLine(:,2),splineLine(:,1),splineLine(:,3),'g','LineWidth',8);


for t=1:length(handles.cartesian)
    im = handles.cartesian(t).Images;
    dim1 = size(im,1);
    dim2 = size(im,2);
    [x,y,z] = meshgrid(1:dim1,1:dim2,1:2);
    rot = handles.cartesian(t).RotationMatrix;
    row1 = nonzeros(rot(1,:));
    row2 = nonzeros(rot(2,:));
    row3 = nonzeros(rot(3,:));
    X = x.*row1(1) + row1(2);
    Y = y.*row2(1) + row2(2);
    Z = z.*row3(1) + row3(2);
    xslice = []; 
    yslice = []; 
    zslice = Z(1,1,1);
    I = repmat(im,[1 1 2]);
    Slice(t) = slice(X,Y,Z,I,xslice,yslice,zslice); 
    handles.Slices = Slice;
    shading interp; colormap gray; hold on;
end 

set(handles.ShowPlanesRadio,'Enable','on','Value',1)



% --- Executes on button press in LoadAxialPush.
function LoadAxialPush_Callback(hObject, eventdata, handles)
[anatomicalFile, anatomicalDir] = uigetfile({'*.dcm;*.dicom;','Useable Files (*.dcm,*.dicom)';
   '*.dcm',  'DICOM files (*.dcm)'; ...
   '*.dicom','DICOM-files (*.dicom)'; ...
   '*.*',  'All Files (*.*)'}, 'Select ONE Axial DICOM image');
[~,~,extension] = fileparts(anatomicalFile); %get file extension
dirInfo = dir(fullfile(anatomicalDir,['*' extension]));
handles.axial.Info = dicominfo(fullfile(anatomicalDir,dirInfo(1).name));
for i=1:length(dirInfo) %read all dcm files
    images(:,:,i) = single(dicomread(fullfile(anatomicalDir,dirInfo(i).name)));
end  
handles.axial.Images = rescale(images);
originShift = [handles.axial.Info.ImagePositionPatient;1]; % origin is top left corner of image
xres = handles.axial.Info.PixelSpacing(1);
yres = handles.axial.Info.PixelSpacing(2);
zres = handles.axial.Info.SliceThickness;

% sometimes get extremely small values that should be 0, so round
xVector = round(handles.axial.Info.ImageOrientationPatient(1:3),8); % what direction rows run w/r/to x
yVector = round(handles.axial.Info.ImageOrientationPatient(4:6),8); % what direction the cols run w/r/to y
zVector = [cross(xVector,yVector);0];

xVector = [xVector;0];
yVector = [yVector;0];
handles.axial.RotationMatrix = ...
    [xres*xVector yres*yVector zres*zVector originShift]; % turn these vectors into matrices

set(handles.ImageSlider,'Enable','on'); 
set(handles.MinContrastUpdate,'Enable','on'); 
set(handles.MaxContrastUpdate,'Enable','on'); 
set(handles.UpdateAxial,'String','Axial Data Loaded'); 
updateAnatImages(handles)

guidata(hObject, handles);



% --- Executes on button press in SegAxialPush.
function SegAxialPush_Callback(hObject, eventdata, handles)
set(handles.ImageSlider,'Enable','off'); 
set(handles.MinContrastUpdate,'Enable','off'); 
set(handles.MaxContrastUpdate,'Enable','off'); 
set(handles.WhiteBoxReminder,'Visible','on');
axes(handles.AnatDisplay); %force axes to anatomical plot

rot = handles.axial.RotationMatrix;
images = handles.axial.Images;
minc = str2double(get(handles.MinContrastUpdate,'String'));
maxc = str2double(get(handles.MaxContrastUpdate,'String'));
x = []; y = []; z = []; %allocate for coordinates of drawn points
for i = 1:size(handles.axial.Images,3) %go slice-by-slice through each image
    im = images(:,:,i);
    leftEndpoint = round(size(im,1).*0.9);
    botEndpoint = round(size(im,2).*0.9);
    largestEndpoint = max(leftEndpoint,botEndpoint);
    im(leftEndpoint:end,botEndpoint:end) = 1;
    imshow(im,[minc maxc]);
    [xTemp,yTemp] = getpts(); %draw points on image along aorta
    zTemp = i.*(ones(size(xTemp,1),1)); %add z-coordinates for slice
    x = [x; xTemp]; %append points
    y = [y; yTemp];
    z = [z; zTemp];
end 

dummy = ones(size(x));
within = x<largestEndpoint;
x = x(within);
y = y(within);
z = z(within);
dummy = dummy(within);
points = [x, y, z, dummy]';

for j=1:size(points,2)
    POINTS(:,j) = rot*points(:,j);
end 
points(4,:) = [];
POINTS(4,:) = [];
handles.axial.points = points;
handles.axial.POINTS = POINTS;
axes(handles.CenterlineDisplay);
h = scatter3(POINTS(1,:),POINTS(2,:),POINTS(3,:), 32, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75], ...
    'DisplayName','Axial');
alpha = 0.3;
set(h,'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
legend('Location','southeast');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');

set(handles.ImageSlider,'Enable','on'); 
set(handles.MinContrastUpdate,'Enable','on'); 
set(handles.MaxContrastUpdate,'Enable','on'); 
set(handles.WhiteBoxReminder,'Visible','off');
guidata(hObject, handles);



% --- Executes on button press in LoadCoronalPush.
function LoadCoronalPush_Callback(hObject, eventdata, handles)
handles.CurrView = 'Coronal';
[anatomicalFile, anatomicalDir] = uigetfile({'*.dcm;*.dicom;','Useable Files (*.dcm,*.dicom)';
   '*.dcm',  'DICOM files (*.dcm)'; ...
   '*.dicom','DICOM-files (*.dicom)'; ...
   '*.*',  'All Files (*.*)'}, 'Select ONE Coronal DICOM image');
[~,~,extension] = fileparts(anatomicalFile); %get file extension
dirInfo = dir(fullfile(anatomicalDir,['*' extension]));
handles.coronal.Info = dicominfo(fullfile(anatomicalDir,dirInfo(1).name));
for i=1:length(dirInfo) %read all dcm files
    images(:,:,i) = single(dicomread(fullfile(anatomicalDir,dirInfo(i).name)));
end  
handles.coronal.Images = rescale(images);
originShift = [handles.coronal.Info.ImagePositionPatient;1]; % origin is top left corner of image
xres = handles.coronal.Info.PixelSpacing(1);
yres = handles.coronal.Info.PixelSpacing(2);
zres = handles.coronal.Info.SliceThickness;

% sometimes get extremely small values that should be 0, so round
xVector = round(handles.coronal.Info.ImageOrientationPatient(1:3),8); % what direction rows run w/r/to x
yVector = round(handles.coronal.Info.ImageOrientationPatient(4:6),8); % what direction the cols run w/r/to y
zVector = [cross(xVector,yVector);0];

xVector = [xVector;0];
yVector = [yVector;0];
handles.coronal.RotationMatrix = ...
    [xres*xVector yres*yVector -zres*zVector originShift]; % turn these vectors into matrices

set(handles.ImageSlider,'Enable','on'); 
set(handles.MinContrastUpdate,'Enable','on'); 
set(handles.MaxContrastUpdate,'Enable','on'); 
set(handles.UpdateCoronal,'String','Coronal Data Loaded'); 
updateAnatImages(handles)

guidata(hObject, handles);



% --- Executes on button press in SegCoronalPush.
function SegCoronalPush_Callback(hObject, eventdata, handles)
set(handles.ImageSlider,'Enable','off'); 
set(handles.MinContrastUpdate,'Enable','off'); 
set(handles.MaxContrastUpdate,'Enable','off'); 
set(handles.WhiteBoxReminder,'Visible','on');
axes(handles.AnatDisplay); %force axes to anatomical plot

rot = handles.coronal.RotationMatrix;
images = handles.coronal.Images;
minc = str2double(get(handles.MinContrastUpdate,'String'));
maxc = str2double(get(handles.MaxContrastUpdate,'String'));
x = []; y = []; z = []; %allocate for coordinates of drawn points
for i = 1:size(handles.coronal.Images,3) %go slice-by-slice through each image
    im = images(:,:,i);
    leftEndpoint = round(size(im,1).*0.9);
    botEndpoint = round(size(im,2).*0.9);
    largestEndpoint = max(leftEndpoint,botEndpoint);
    im(leftEndpoint:end,botEndpoint:end) = 1;
    imshow(im,[minc maxc]);
    [xTemp,yTemp] = getpts(); %draw points on image along aorta
    zTemp = i.*(ones(size(xTemp,1),1)); %add z-coordinates for slice
    x = [x; xTemp]; %append points
    y = [y; yTemp];
    z = [z; zTemp];
end 

dummy = ones(size(x));
within = x<largestEndpoint;
x = x(within);
y = y(within);
z = z(within);
dummy = dummy(within);
points = [x, y, z, dummy]';

for j=1:size(points,2)
    POINTS(:,j) = rot*points(:,j);
end 
points(4,:) = [];
POINTS(4,:) = [];
handles.coronal.points = points;
handles.coronal.POINTS = POINTS;
axes(handles.CenterlineDisplay); hold on;
h = scatter3(POINTS(1,:),POINTS(2,:),POINTS(3,:), 32, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.75 .75 0], ...
    'DisplayName','Coronal');
alpha = 0.3;
set(h,'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
legend('Location','southeast');

set(handles.ImageSlider,'Enable','on'); 
set(handles.MinContrastUpdate,'Enable','on'); 
set(handles.MaxContrastUpdate,'Enable','on'); 
set(handles.WhiteBoxReminder,'Visible','off');
guidata(hObject, handles);



% --- Executes on button press in LoadSagittalPush.
function LoadSagittalPush_Callback(hObject, eventdata, handles)
handles.CurrView = 'Sagittal';
[anatomicalFile, anatomicalDir] = uigetfile({'*.dcm;*.dicom;','Useable Files (*.dcm,*.dicom)';
   '*.dcm',  'DICOM files (*.dcm)'; ...
   '*.dicom','DICOM-files (*.dicom)'; ...
   '*.*',  'All Files (*.*)'}, 'Select ONE Sagittal DICOM image');
[~,~,extension] = fileparts(anatomicalFile); %get file extension
dirInfo = dir(fullfile(anatomicalDir,['*' extension]));
handles.sagittal.Info = dicominfo(fullfile(anatomicalDir,dirInfo(1).name));
for i=1:length(dirInfo) %read all dcm files
    images(:,:,i) = single(dicomread(fullfile(anatomicalDir,dirInfo(i).name)));
end  
handles.sagittal.Images = rescale(images);
originShift = [handles.sagittal.Info.ImagePositionPatient;1]; % origin is top left corner of image
xres = handles.sagittal.Info.PixelSpacing(1);
yres = handles.sagittal.Info.PixelSpacing(2);
zres = handles.sagittal.Info.SliceThickness;

% sometimes get extremely small values that should be 0, so round
xVector = round(handles.sagittal.Info.ImageOrientationPatient(1:3),8); % what direction rows run w/r/to x
yVector = round(handles.sagittal.Info.ImageOrientationPatient(4:6),8); % what direction the cols run w/r/to y
zVector = [cross(xVector,yVector);0];

xVector = [xVector;0];
yVector = [yVector;0];
handles.sagittal.RotationMatrix = ...
    [xres*xVector yres*yVector zres*zVector originShift]; % turn these vectors into matrices

set(handles.ImageSlider,'Enable','on'); 
set(handles.MinContrastUpdate,'Enable','on'); 
set(handles.MaxContrastUpdate,'Enable','on'); 
set(handles.UpdateSagittal,'String','Sagittal Data Loaded'); 
updateAnatImages(handles)

guidata(hObject, handles);



% --- Executes on button press in SegSagittalPush.
function SegSagittalPush_Callback(hObject, eventdata, handles)
set(handles.ImageSlider,'Enable','off'); 
set(handles.MinContrastUpdate,'Enable','off'); 
set(handles.MaxContrastUpdate,'Enable','off'); 
set(handles.WhiteBoxReminder,'Visible','on');
axes(handles.AnatDisplay); %force axes to anatomical plot

rot = handles.sagittal.RotationMatrix;
images = handles.sagittal.Images;
minc = str2double(get(handles.MinContrastUpdate,'String'));
maxc = str2double(get(handles.MaxContrastUpdate,'String'));
x = []; y = []; z = []; %allocate for coordinates of drawn points
for i = 1:size(handles.sagittal.Images,3) %go slice-by-slice through each image
    im = images(:,:,i);
    leftEndpoint = round(size(im,1).*0.9);
    botEndpoint = round(size(im,2).*0.9);
    largestEndpoint = max(leftEndpoint,botEndpoint);
    im(leftEndpoint:end,botEndpoint:end) = 1;
    imshow(im,[minc maxc]);
    [xTemp,yTemp] = getpts(); %draw points on image along aorta
    zTemp = i.*(ones(size(xTemp,1),1)); %add z-coordinates for slice
    x = [x; xTemp]; %append points
    y = [y; yTemp];
    z = [z; zTemp];
end 

dummy = ones(size(x));
within = x<largestEndpoint;
x = x(within);
y = y(within);
z = z(within);
dummy = dummy(within);
points = [x, y, z, dummy]';

for j=1:size(points,2)
    POINTS(:,j) = rot*points(:,j);
end 
points(4,:) = [];
POINTS(4,:) = [];
handles.sagittal.points = points;
handles.sagittal.POINTS = POINTS;
axes(handles.CenterlineDisplay); hold on;
h = scatter3(POINTS(1,:),POINTS(2,:),POINTS(3,:), 32, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.75 0 .75], ...
    'DisplayName','Sagittal');
alpha = 0.3;
set(h,'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
legend('Location','southeast');

set(handles.ImageSlider,'Enable','on'); 
set(handles.MinContrastUpdate,'Enable','on'); 
set(handles.MaxContrastUpdate,'Enable','on'); 
set(handles.WhiteBoxReminder,'Visible','off');
guidata(hObject, handles);



% --- Executes on button press in Seg2DCartPush.
function Seg2DCartPush_Callback(hObject, eventdata, handles)
handles.CurrView = '2DCartesian';
cartIter = handles.CurrCartesian;
[anatomicalFile, anatomicalDir] = uigetfile({'*.dcm;*.dicom;','Useable Files (*.dcm,*.dicom)';
   '*.dcm',  'DICOM files (*.dcm)'; ...
   '*.dicom','DICOM-files (*.dicom)'; ...
   '*.*',  'All Files (*.*)'}, 'Select ONE 2DPC Cartesian DICOM image');
[~,~,extension] = fileparts(anatomicalFile); %get file extension
dirInfo = dir(fullfile(anatomicalDir,['*' extension]));
handles.cartesian(cartIter).Info = dicominfo(fullfile(anatomicalDir,dirInfo(1).name));
for i=1:length(dirInfo) %read all dcm files
    images(:,:,i) = single(dicomread(fullfile(anatomicalDir,dirInfo(i).name)));
end  
mag = images(:,:,(size(images,3)/2)+1:end);
MAG = mean(mag,3);

handles.cartesian(cartIter).Images = rescale(MAG);
originShift = [handles.cartesian(cartIter).Info.ImagePositionPatient;1]; % origin is top left corner of image
xres = handles.cartesian(cartIter).Info.PixelSpacing(1);
yres = handles.cartesian(cartIter).Info.PixelSpacing(2);
zres = handles.cartesian(cartIter).Info.SliceThickness;

% sometimes get extremely small values that should be 0, so round
xVector = round(handles.cartesian(cartIter).Info.ImageOrientationPatient(1:3),8); % what direction rows run w/r/to x
yVector = round(handles.cartesian(cartIter).Info.ImageOrientationPatient(4:6),8); % what direction the cols run w/r/to y
zVector = [cross(xVector,yVector);0];

xVector = [xVector;0];
yVector = [yVector;0];
handles.cartesian(cartIter).RotationMatrix = ...
[xres*xVector yres*yVector zres*zVector originShift]; % turn vectors into matrices


axes(handles.AnatDisplay); %force axes to anatomical plot

rot = handles.cartesian(cartIter).RotationMatrix;
minc = str2double(get(handles.MinContrastUpdate,'String'));
maxc = str2double(get(handles.MaxContrastUpdate,'String'));

imshow(handles.cartesian(cartIter).Images,[minc maxc]);
[x,y] = getpts(); %draw points on image along aorta
z = (ones(size(x,1),1)); %add z-coordinates for slice
dummy = ones(size(x));
points = [x, y, z, dummy]';

for j=1:size(points,2)
    POINTS(:,j) = rot*points(:,j);
end 
points(4,:) = [];
POINTS(4,:) = [];
handles.cartesian(cartIter).points = points;
handles.cartesian(cartIter).POINTS = POINTS;
axes(handles.CenterlineDisplay); hold on;
scatter3(POINTS(1,:),POINTS(2,:),POINTS(3,:),'k*', ...
    'LineWidth',12, ...
    'DisplayName','Cart-2DPC');
legend('Location','southeast','AutoUpdate','off');

handles.CurrCartesian = cartIter + 1;
guidata(hObject, handles);



% --- Executes on button press in MinContrastUpdate.
function MinContrastUpdate_Callback(hObject, eventdata, handles)
updateAnatImages(handles)

% --- Executes during object creation, after setting all properties.
function MinContrastUpdate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in MaxContrastUpdate.
function MaxContrastUpdate_Callback(hObject, eventdata, handles)
updateAnatImages(handles)

% --- Executes during object creation, after setting all properties.
function MaxContrastUpdate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in CreateCenterlinePush.
function CreateCenterlinePush_Callback(hObject, eventdata, handles)
PC = [];
for p=1:length(handles.cartesian)
    PC = [PC handles.cartesian(p).POINTS];
end 
figure; 
scatter(handles.coronal.POINTS(1,:),handles.coronal.POINTS(3,:),'b'); hold on;
scatter(handles.axial.POINTS(1,:),handles.axial.POINTS(3,:),'b'); hold on;
scatter(PC(1,:),PC(3,:),'filled','g');
title('Coronal Projection');
cor = drawpolyline;
corP = cor.Position;
handles.TraceCoronal = gcf;

figure; scatter(handles.sagittal.POINTS(2,:),handles.sagittal.POINTS(3,:),'b');
hold on; scatter(handles.axial.POINTS(2,:),handles.axial.POINTS(3,:),'b');
hold on; scatter(PC(2,:),PC(3,:),'filled','g');
title('Sagittal Projection');
sag = drawpolyline;
sagP = sag.Position;
handles.TraceSagittal = gcf;

splinePositions1 = interppolygon(sagP,150); %interpolate sagittal polyline (150 points)
splinePositions2 = interppolygon(corP,150); %interpolate coronal polyline (150 points)
splineLine(:,1) = splinePositions1(:,1);
splineLine(:,2) = splinePositions2(:,1);
splineLine(:,3) = splinePositions1(:,2);
axes(handles.CenterlineDisplay);
plot3(splineLine(:,2),splineLine(:,1),splineLine(:,3),'g','LineWidth',8);

for t=1:length(handles.cartesian)
    im = handles.cartesian(t).Images;
    dim1 = size(im,1);
    dim2 = size(im,2);
    [x,y,z] = meshgrid(1:dim1,1:dim2,1:2);
    rot = handles.cartesian(t).RotationMatrix;
    row1 = nonzeros(rot(1,:));
    row2 = nonzeros(rot(2,:));
    row3 = nonzeros(rot(3,:));
    X = x.*row1(1) + row1(2);
    Y = y.*row2(1) + row2(2);
    Z = z.*row3(1) + row3(2);
    xslice = []; 
    yslice = []; 
    zslice = Z(1,1,1);
    I = repmat(im,[1 1 2]);
    Slice(t) = slice(X,Y,Z,I,xslice,yslice,zslice); 
    handles.Slices = Slice;
    shading interp; colormap gray; hold on;
end 

set(handles.ShowPlanesRadio,'Enable','on','Value',1)
handles.Centerline = splineLine;
guidata(hObject, handles);



% --- Executes on button press in SaveCenterlinePush.
function SaveCenterlinePush_Callback(hObject, eventdata, handles)
centerline = handles.Centerline;
PC = [];
for p=1:length(handles.cartesian)
    PC = [PC handles.cartesian(p).POINTS];
end 

for r=1:size(PC,2)
    diff = centerline-PC(:,r)';
    dist = vecnorm(diff');
    [~,idx] = min(dist);
    IDX(r) = idx;
end 

distances = zeros(1,size(centerline,1));
for d=2:size(centerline,1)
    distances(d) = sqrt( sum( (centerline(d-1,:)-centerline(d,:)).^2 ) );
end 

for i=1:(length(IDX)-1)
    PlaneDistances(i) = sum( distances(IDX(i):(IDX(i+1)-1)) );
end 

anatCLdataset.Axial = handles.axial;
anatCLdataset.Coronal = handles.coronal;
anatCLdataset.Sagittal = handles.sagittal;
anatCLdataset.Cartesian2DPC = handles.cartesian;
anatCLdataset.Centerline = centerline;
anatCLdataset.Distances = distances;
anatCLdataset.ROIindices = IDX;
anatCLdataset.PlaneDistances = PlaneDistances;
directory = uigetdir(); %saving location
save([directory filesep 'anatCLdataset.mat'],'anatCLdataset');

frame = getframe(handles.CenterlineDisplay);
imwrite(frame2im(frame),[directory filesep 'CenterlineData' filesep 'Centerline3D.png']);
frame = getframe(handles.TraceCoronal);
imwrite(frame2im(frame),[directory filesep 'CenterlineData' filesep 'CoronalTrace.png']);
frame = getframe(handles.TraceSagittal);
imwrite(frame2im(frame),[directory filesep 'CenterlineData' filesep 'SagittalTrace.png']);

set(handles.SaveText,'String','Centerline Saved!'); 



% --- Executes on slider movement.
function ImageSlider_Callback(hObject, eventdata, handles)
updateAnatImages(handles)

% --- Executes during object creation, after setting all properties.
function ImageSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in ShowPlanesRadio.
function ShowPlanesRadio_Callback(hObject, eventdata, handles)
if get(handles.ShowPlanesRadio,'Value')==1
    Slices = handles.Slices;
    for t=1:length(Slices)
        axes(handles.CenterlineDisplay);
        set(Slices(t),'visible','on');
    end 
else
    Slices = handles.Slices;
    for t=1:length(Slices)
        axes(handles.CenterlineDisplay);
        set(Slices(t),'visible','off');
    end 
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Update Images in ANATOMICAL PLOT
function updateAnatImages(handles)
    axes(handles.AnatDisplay); %force axes to anatomical plot
    
    switch handles.CurrView
        case 'Axial'
            images = handles.axial.Images;
        case 'Coronal'
            images = handles.coronal.Images;
        case 'Sagittal'
            images = handles.sagittal.Images;
        case '2DCartesian'
            images = handles.cartesian.Images;
        case '2DRadial'
            images = handles.radial.Images;
    end 
    
    minc = str2double(get(handles.MinContrastUpdate,'String'));
    maxc = str2double(get(handles.MaxContrastUpdate,'String'));

    if ndims(images)<3 %if we only have one image (should rarely happen)
        maxSize = max(size(images,1),size(images,2));
        steps = [1 maxSize]; %set our slider to as wide as possible 
        set(handles.ImageSlider,'SliderStep', steps);
        anatSlice = images;
    else
        dim3size = size(images,3);
        steps = [1/(dim3size-1) 10/(dim3size-1)];
        set(handles.ImageSlider,'SliderStep', steps);
        sliceNum = 1+round( get(handles.ImageSlider,'Value').*(dim3size-1) );  %get slice number from slider
        anatSlice = images(:,:,sliceNum); %get slice from images
    end 
    imshow(anatSlice,[minc maxc]); %show image
return
    


% --- Turn PolyLine into SplineLine    
function Y = interppolygon(X,N)
    if nargin < 2 || N < 2
        N = 2; %if only one arg or too small N, just assume 2
    end
    nDim = size(X,2); %should be 2
    dx = 0;
    
    for dim = 1:nDim
        dx = dx + diff(X(:,dim)).^2 ; %get sum of squares in each dim
    end
    
    lengthBetweenPoints = sqrt(dx); %now get distance
    lengthLine = sum(lengthBetweenPoints);
    origMetric = [0; cumsum(lengthBetweenPoints/lengthLine)];
    
    interpMetric = (0:(1/(N-1)):1)';
    Y = interp1(origMetric,X,interpMetric,'makima'); %makima seems to work well
    %Y = csaps([0 times times(end)+times(1)],[0 meanROI 0],0.0001,timesInterp);
return



function v = load_dat(name, matrix)
[fid,errmsg]= fopen(name,'r');
if fid < 0  % If name does not exist in directory
    disp(['Error Opening Data : ',errmsg]);
end

% Reads in as short, reshapes by image resolution (i.e. 320x320x320)
v = reshape(fread(fid,'short=>single'),matrix);
fclose(fid);

return



% --- Executes on button press in Seg2DRadialPush.
% function Seg2DRadialPush_Callback(hObject, eventdata, handles)
% handles.CurrView = '2DRadial';
% radIter = handles.CurrRadial;
% 
% [~, anatomicalDir] = uigetfile({'*.dat','Useable Files (*.dat)';
%    '*.dat',  'DAT files (*.dat)'; ...
%    '*.*',  'All Files (*.*)'}, 'Select ONE 2DPC Cartesian DAT dataset');
% 
% fid = fopen([anatomicalDir 'pcvipr_header.txt'], 'r'); %open header
% dataArray = textscan(fid,'%s%s%[^\n\r]','Delimiter',' ', ...
%     'MultipleDelimsAsOne',true,'ReturnOnError',false); %parse header info
% fclose(fid);
% dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:),'UniformOutput',false);
% pcviprHeader = cell2struct(dataArray{1,2}(:),dataArray{1,1}(:),1); %turn to structure
% handles.radial(radIter).Info = pcviprHeader; %add pcvipr header to handles
% dimx = pcviprHeader.matrixx; %resolution in x
% dimy = pcviprHeader.matrixy; %resolution in y
% 
% MAG = load_dat(fullfile(anatomicalDir,'MAG.dat'),[dimx dimy]); %Average magnitude
% MAG = flipud(MAG);
% handles.radial(radIter).Images = rescale(MAG);
% 
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
% handles.radial(radIter).RotationMatrix = ...
%     [xVector yVector zVector originShift];
% 
% 
% axes(handles.AnatDisplay); %force axes to anatomical plot
% 
% rot = handles.radial(radIter).RotationMatrix;
% minc = str2double(get(handles.MinContrastUpdate,'String'));
% maxc = str2double(get(handles.MaxContrastUpdate,'String'));
% 
% imshow(handles.radial(radIter).Images,[minc maxc]);
% [x,y] = getpts(); %draw points on image along aorta
% z = (ones(size(x,1),1)); %add z-coordinates for slice
% dummy = ones(size(x));
% points = [x, y, z, dummy]';
% 
% for j=1:size(points,2)
%     POINTS(:,j) = rot*points(:,j);
% end 
% points(4,:) = [];
% POINTS(4,:) = [];
% handles.cartesian(radIter).points = points;
% handles.cartesian(radIter).POINTS = POINTS;
% axes(handles.CenterlineDisplay); hold on;
% scatter3(POINTS(1,:),POINTS(2,:),POINTS(3,:),'r*', ...
%     'LineWidth',12, ...
%     'DisplayName','Rad-2DPC');
% legend('Location','southeast');
% 
% handles.CurrRadial = radIter + 1;
% guidata(hObject, handles);

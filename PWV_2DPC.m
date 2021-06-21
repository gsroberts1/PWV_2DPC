function varargout = PWV_2DPC(varargin)
    % PWV_2DPC MATLAB code for PWV_2DPC.fig
    %      PWV_2DPC, by itself, creates a new PWV_2DPC or raises the existing
    %      singleton*.
    %
    %      H = PWV_2DPC returns the handle to a new PWV_2DPC or the handle to
    %      the existing singleton*.
    %
    %      PWV_2DPC('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in PWV_2DPC.M with the given input arguments.
    %
    %      PWV_2DPC('Property','Value',...) creates a new PWV_2DPC or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before PWV_2DPC_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to PWV_2DPC_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help PWV_2DPC

    % Last Modified by GUIDE v2.5 03-Jun-2021 15:17:10
    % Developed by Grant S Roberts, University of Wisconsin-Madison, 2019
    
    
    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @PWV_2DPC_OpeningFcn, ...
                       'gui_OutputFcn',  @PWV_2DPC_OutputFcn, ...
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


function PWV_2DPC_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to PWV_2DPC (see VARARGIN)

    % Choose default command line output for PWV_2DPC
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % Get anatomical, pc, and bSSFP data from LoadPWV GUI
    handles.pcDatasets(1).ROI = []; %initialize ROI field
    handles.pcDatasets(1).ROIdataUninterp = []; %throgh-plane flow curve data
    handles.pcDatasets(1).ROIdataGaussian = []; %flow smoothed with Gauss.
    
    handles.global.interpType = 'Gaussian'; %interpolation type (i.e. Gaussian)
    handles.global.startAnalyzing = 0; %flag to begin PWV calculations
    handles.global.totalROIs = 0;
    handles.global.pgShift = 0;
    handles.global.pcIter = 1;

    set(handles.load2DPCbutton,'Enable','off');
    set(handles.pcPlanePopup,'Enable','off');
    set(handles.pcDatasetPopup,'Enable','off');
    set(handles.drawROIbutton,'Enable','off');
    set(handles.loadROIbutton,'Enable','off');
    set(handles.pcSlider,'Enable','off');
    set(handles.interpolatePopup,'String',{'Gaussian','None'}); %set all possible interpolation types
    set(handles.interpolatePopup,'Enable','off'); %initialize radios and buttons
    set(handles.errorBarRadio,'Enable','off');
    set(handles.pgShiftRadio,'Enable','off');
    set(handles.ttpointRadio,'Value',1); 
    set(handles.ttpointRadio,'Enable','off');
    set(handles.ttuRadio,'Value',1);
    set(handles.ttuRadio,'Enable','off');
    set(handles.ttfRadio,'Value',1);
    set(handles.ttfRadio,'Enable','off');
    set(handles.xcorrRadio,'Value',1);
    set(handles.xcorrRadio,'Enable','off');
    set(handles.exportAnalysisButton,'Enable','off');
    
    guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = PWV_2DPC_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;


    
%%%%%%%%%%%% LOAD 2DPC PLANE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PLANE PLOT - CREATE FUNCTION
function pcPlanePlot_CreateFcn(hObject, eventdata, handles)


% --- LOAD CENTERLINE DATA - CALLBACK
function loadCLpush_Callback(hObject, eventdata, handles)
    [clFile, clDir] = uigetfile({'*.mat;','Useable Files (*.mat)';
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select the centerline dataset (anatCLdataset.mat)');
    load([clDir clFile]);
    handles.centerline = anatCLdataset;
    
    cd(clDir);
    handles.global.homeDir = clDir;
    set(handles.load2DPCbutton,'Enable','on');
    guidata(hObject, handles);
    
    
% --- LOAD 2DPC DATASETS - CALLBACK
function load2DPCbutton_Callback(hObject, eventdata, handles)
    [pcFile, pcDir] = uigetfile({'*.dcm;*.dat;*.mat','Useable Files (*.dcm,*.dat,*.mat)';
       '*.dcm',  'DICOM files (*.dcm)'; ...
       '*.dat',  'DAT-files (*.dat)'; ...
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select ONE 2DPC file in the dataset');
    pcIter = handles.global.pcIter;
    [~,~,extension] = fileparts(pcFile);
    dirInfo = dir(fullfile(pcDir,['*' extension]));
    if isequal(extension,'.dcm') %if our extension is a dicom file
        handles.pcDatasets(pcIter).Info = dicominfo(fullfile(pcDir,dirInfo(1).name)); %get dicom metadata (from 1st dicom)
        for i=1:length(dirInfo)
            hold(:,:,i) = single(dicomread(fullfile(pcDir,dirInfo(i).name))); %read dicoms and cast to single
        end  
        mag = hold(:,:,floor(length(dirInfo)/2)+1:end); %magnitude is last half
        v = hold(:,:,1:floor(length(dirInfo)/2)); %velocity is first half of images
        MAG = mean(mag,3); %time-averaged magnitude
        VMEAN = mean(v,3); %time-averaged velocity
        CD = MAG.*sin( pi/2*abs(VMEAN)/max(VMEAN(:)) );
        
        handles.pcDatasets(pcIter).Images.MAG = MAG;
        handles.pcDatasets(pcIter).Images.CD = CD;
        handles.pcDatasets(pcIter).Images.V = VMEAN;
        handles.pcDatasets(pcIter).Images.mag = mag;
        handles.pcDatasets(pcIter).Images.v = v; 
        handles.pcDatasets(pcIter).Names = ['Plane' num2str(pcIter)];
    elseif isequal(extension,'.dat')
        fid = fopen([pcDir filesep 'pcvipr_header.txt'], 'r'); %open header
        dataArray = textscan(fid,'%s%s%[^\n\r]','Delimiter',' ', ...
            'MultipleDelimsAsOne',true,'ReturnOnError',false); %parse header info
        fclose(fid);
        dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:),'UniformOutput',false);
        pcviprHeader = cell2struct(dataArray{1,2}(:),dataArray{1,1}(:),1); %turn to structure
        handles.pcDatasets(pcIter).Info = pcviprHeader; %add pcvipr header to handles
        resx = pcviprHeader.matrixx; %resolution in x
        resy = pcviprHeader.matrixy; %resolution in y
        nframes = pcviprHeader.frames; %number of cardiac frames
        MAG = load_dat(fullfile(pcDir,'MAG.dat'),[resx resy]); %Average magnitude
        CD = load_dat(fullfile(pcDir,'CD.dat'),[resx resy]); %Average complex difference
        VMEAN = load_dat(fullfile(pcDir,'comp_vd_3.dat'),[resx resy]); %Average velocity

        % Initialize data time-resolved data arrays
        mag = zeros(resx,resy,nframes); %Time-resolved magnitude
        cd = zeros(resx,resy,nframes); %Time-resolved complex difference
        v = zeros(resx,resy,nframes); %Time-resolved velocity 
        for j = 1:nframes  %velocity is placed in v3 for 2D (through-plane)
            mag(:,:,j) = load_dat(fullfile(pcDir,[filesep 'ph_' num2str(j-1,'%03i') '_mag.dat']),[resx resy]);
            cd(:,:,j) = load_dat(fullfile(pcDir,[filesep 'ph_' num2str(j-1,'%03i') '_cd.dat']),[resx resy]);
            v(:,:,j) = load_dat(fullfile(pcDir,[filesep 'ph_' num2str(j-1,'%03i') '_vd_3.dat']),[resx resy]);
        end
        
        handles.pcDatasets(pcIter).Images.MAG = flipud(MAG);
        handles.pcDatasets(pcIter).Images.CD = flipud(CD);
        handles.pcDatasets(pcIter).Images.V = flipud(VMEAN);
        handles.pcDatasets(pcIter).Images.mag = flipud(mag);
        handles.pcDatasets(pcIter).Images.cd = flipud(cd);
        handles.pcDatasets(pcIter).Images.v = flipud(v); 
        handles.pcDatasets(pcIter).Names = ['Plane' num2str(pcIter)];
    else %if a single matlab file (with all images)
        hold = load([pcDir pcFile]);
        planeName = fieldnames(hold);
        planeName = planeName{1};
        handles.pcDatasets(pcIter).Info = hold.(planeName).Info;
        images = hold.(planeName).Images;
        images = single(images);
        
        handles.pcDatasets(pcIter).Images.MAG = mean(images(:,:,:,1),3);
        handles.pcDatasets(pcIter).Images.CD = mean(images(:,:,:,2),3);
        handles.pcDatasets(pcIter).Images.V = mean(images(:,:,:,3),3);
        handles.pcDatasets(pcIter).Images.mag = images(:,:,:,1);
        handles.pcDatasets(pcIter).Images.cd = images(:,:,:,2);
        handles.pcDatasets(pcIter).Images.v = images(:,:,:,3);
        handles.pcDatasets(pcIter).Names = planeName;
    end
    
    handles.global.pcIter = handles.global.pcIter + 1;
    
    set(handles.pcPlanePopup,'Enable','on');
    set(handles.pcDatasetPopup,'Enable','on');
    set(handles.drawROIbutton,'Enable','on');
    set(handles.loadROIbutton,'Enable','on');
    set(handles.pcSlider,'Enable','on');
    set(handles.pcPlanePopup,'String',{handles.pcDatasets.Names}); %list of all planes (AAo, AbdAo, etc.)
    set(handles.pcDatasetPopup,'String',fieldnames(handles.pcDatasets(pcIter).Images)); %list of all datasets (CD, MAG, v, etc.)
    guidata(hObject, handles);
    updatePCImages(handles);
    
    
% --- PLANE DROPDOWN - CALLBACK
function pcPlanePopup_Callback(hObject, eventdata, handles)   
    updatePCImages(handles); %update images on PC plot anytime we click on a new plane

% --- PLANE DROPDOWN - CREATE FUNCTION
function pcPlanePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- DATASET DROPDOWN - CALLBACK
function pcDatasetPopup_Callback(hObject, eventdata, handles)
    updatePCImages(handles); %update images on PC plot anytime we click on a new dataset

% --- DATASET DROPDOWN - CREATE FUNCTION
function pcDatasetPopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- DRAWROI BUTTON - CALLBACK
function drawROIbutton_Callback(hObject, eventdata, handles)
    axes(handles.pcPlanePlot); %make sure we do this on top PC plot 

    set(handles.loadCLpush,'Enable','off');
    set(handles.load2DPCbutton,'Enable','off');
    set(handles.pcPlanePopup,'Enable','off'); %make it so we can't select another plane
    set(handles.pcDatasetPopup,'Enable','off'); %make it so we can't select another plane
    set(handles.drawROIbutton,'Enable','off'); %make it so we can't draw a second ROI
    
    planeNum = get(handles.pcPlanePopup,'Value'); %get current plane of interest (eg AAo)
    mydlg = warndlg('Press enter when the ROI is set'); %open dialog warning 
    waitfor(mydlg); %MUST PRESS ENTER TO PROCEED
    circle = drawcircle('FaceAlpha',0.1,'Color','g','LineWidth',1,'Deletable',0); %draw circle on PC image
    while true
        w = waitforbuttonpress; %wait for enter push ...
        switch w 
            case 1 % if it was a keyboard press.
            key = get(gcf,'currentcharacter'); %get key that was pressed
                switch key
                    case 27 % escape key
                        set(handles.loadCLpush,'Enable','on');
                        set(handles.load2DPCbutton,'Enable','on');
                        set(handles.pcPlanePopup,'Enable','on'); %make it so we can't select another plane
                        set(handles.pcDatasetPopup,'Enable','on'); %make it so we can't select another plane
                        set(handles.drawROIbutton,'Enable','on');
                        broke = 1;
                        break % break out of the while loop
                    case 13 % 13 is the enter/return key 
                        circle.InteractionsAllowed = 'none'; %freeze circle
                        broke = 0;
                        break
                    otherwise 
                        %wait for a different command
                end
       end
    end
    
    if ~broke
        handles.pcDatasets(planeNum).ROI = circle; %temporarily hold circle (deleted once gone)
    else 
        handles.pcDatasets(planeNum).ROI = [];
    end 

    guidata(hObject,handles);
    updatePCImages(handles); 

    
% --- LOAD ROI BUTTON - CALLBACK
function loadROIbutton_Callback(hObject, eventdata, handles)
    planeNum = get(handles.pcPlanePopup,'Value'); %get current plane (eg AAo)
    handles.global.totalROIs = handles.global.totalROIs + 1; %add 1 to total ROI count

    v = handles.pcDatasets(planeNum).Images.v; %grab time-resolved velocity
    circle = handles.pcDatasets(planeNum).ROI; %pull circle data from handles
    radius = circle.Radius; %get radius of circle
    center = round(circle.Center); %get center coordinates

    [X,Y] = ndgrid(1:size(v,1),1:size(v,2));
    X = X-center(2); %shift coordinate grid
    Y = Y-center(1);
    roiMask = sqrt(X.^2 + Y.^2)<=radius; %anything outside radius is ignored

    %%% Create Linear Interpolated Data
    if isfield(handles.pcDatasets(planeNum).Info,'matrixx') %if radial data (pcvipr recon)
        matrixx = handles.pcDatasets(planeNum).Info.matrixx; %matrix size in x dimension
        fovx = handles.pcDatasets(planeNum).Info.fovx;  %field of view (mm)
        xres = fovx/matrixx; %resolution (mm). ASSUMED TO BE SAME IN Y DIMENSION
%         frames = handles.pcDatasets(planeNum).Info.frames;
%         timeres = handles.pcDatasets(planeNum).Info.timeres; %temporal resolution (ms)
        frames = handles.pcDatasets(1).Info.frames;
        timeres = handles.pcDatasets(1).Info.timeres;
    else 
        xres = handles.pcDatasets(planeNum).Info.PixelSpacing(1); %resolution (mm) ASSUMED SAME IN Y DIM
%         rrInterval = handles.pcDatasets(planeNum).Info.NominalInterval; %average RR interval (ms)
%         frames = handles.pcDatasets(planeNum).Info.CardiacNumberOfImages; %number of cardiac frames
%         timeres = rrInterval/frames; %temporal resolution (ms)
        rrInterval = handles.pcDatasets(1).Info.NominalInterval; %avg RR int. (ms) TAKE FIRST TO MATCH TEMP RES
        frames = handles.pcDatasets(planeNum).Info.CardiacNumberOfImages; %number of cardiac frames
        timeres = rrInterval/frames; %temporal resolution (ms)
    end 

    area = sum(roiMask(:))*(xres)^2; %ROI area (mm^2)
    for i=1:frames
        vTemp = v(:,:,i); %through-plane velocity in frame i
        roiDataRaw(:,i) = double(vTemp(roiMask)); %indexed velocities within mask
        meanROI(i) = mean(roiDataRaw(:,i)); %mean velocity in frame i (mm/s)
        stdvROI(i) = std(double(vTemp(roiMask))); %stdv of velocity in frame i (mm/s)
        flowROI(i) = area.*meanROI(i).*0.001; %flow in frame i (mm^3/s = mL/s)
    end 

    times = double(timeres.*(0:(frames-1))); %original times
%     tq = 0:0.1:size(v,3)-1; %interpolate time dimension
    sampleDensity = 1000;
    tq = linspace(0,frames-1,sampleDensity); %interpolate time dimension
    timesInterp = double(timeres.*tq); %interpolated times

    %Linear interpolation (to get more points on flow curve)
    meanROIfit = interp1(times,meanROI,timesInterp,'linear');
    stdvROIfit = interp1(times,stdvROI,timesInterp,'linear');
    flowROIfit = interp1(times,flowROI,timesInterp,'linear');

    %Add data to roiStatistics structure
    roiInfo.radius = radius; 
    roiInfo.center = center;
    roiInfo.roiMask = roiMask; 
    roiInfo.roiDataRaw = roiDataRaw;
    roiInfo.Name = ''; %will get changed below
    roiInfo.ROInumber = handles.global.totalROIs;
    
    roiStatisticsUninterp.times = timesInterp;
    roiStatisticsUninterp.meanROI = meanROIfit; 
    roiStatisticsUninterp.stdvROI = stdvROIfit;
    roiStatisticsUninterp.flowROI = flowROIfit; 


    %%% Create Interpolated Curve with Gaussian Smoothing           
    meanROIfit = interp1(times,smoothdata(meanROI,'gaussian',3),timesInterp,'cubic');
    stdvROIfit = interp1(times,smoothdata(stdvROI,'gaussian',3),timesInterp,'cubic');
    flowROIfit = interp1(times,smoothdata(flowROI,'gaussian',3),timesInterp,'cubic');

    roiStatisticsGaussian.times = timesInterp;
    roiStatisticsGaussian.meanROI = meanROIfit; 
    roiStatisticsGaussian.stdvROI = stdvROIfit; 
    roiStatisticsGaussian.flowROI = flowROIfit; 

    %%% Save all data into handles
    if isstruct(handles.pcDatasets(planeNum).ROIdataUninterp)
        handles.pcDatasets(planeNum).ROIinfo(end+1) = roiInfo;
        handles.pcDatasets(planeNum).ROIdataUninterp(end+1) = roiStatisticsUninterp;
        handles.pcDatasets(planeNum).ROIdataGaussian(end+1) = roiStatisticsGaussian;
    else
        handles.pcDatasets(planeNum).ROIinfo = roiInfo;
        handles.pcDatasets(planeNum).ROIdataUninterp = roiStatisticsUninterp;
        handles.pcDatasets(planeNum).ROIdataGaussian = roiStatisticsGaussian;
    end 

    dataDir = handles.global.homeDir; %directory in which plane data is located
    if dataDir(end)=='\' || dataDir(end)=='/' %kill the slash if it exists
        dataDir(end) = [];
    end 
    if ~exist([dataDir filesep 'ROIimages'],'dir') %if the proposed directory doesn't exist
        mkdir([dataDir filesep 'ROIimages']); %make it
        cd([dataDir filesep 'ROIimages']); %move into it
        frame = getframe(handles.pcPlanePlot); %get a snapshot of the PC plane plot with ROI
        image = frame2im(frame); %make into image
        imwrite(image,[handles.pcDatasets(planeNum).Names '.png']) %write it out as PNG
    else
        cd([dataDir filesep 'ROIimages']); %if ROIimages already exists, move into it
        frame = getframe(handles.pcPlanePlot);
        image = frame2im(frame);
        imwrite(image,[handles.pcDatasets(planeNum).Names '.png'])
    end 
    cd(handles.global.homeDir); %lets go back home
    
    %%% Label each ROI w/ names (helpful because there may be 2 ROIs/plane)
    for i=1:numel(handles.pcDatasets)
        if isstruct(handles.pcDatasets(i).ROIinfo) %if we've made ROI data for this dataset
            if length(handles.pcDatasets(i).ROIinfo)==1
                handles.pcDatasets(i).ROIinfo.Name = handles.pcDatasets(planeNum).Names; %name ROI
            else
                for j=1:length(handles.pcDatasets(i).ROIinfo)
                    name = handles.pcDatasets(i).Names; %get plane name
                    planeName = [name ' ROI ' num2str(j)]; %needed if more than one ROI/plane
                    handles.pcDatasets(i).ROIinfo(j).Name = planeName; %name ROI
                end 
            end 
        end 
    end 

    set(handles.interpolatePopup,'Enable','on'); %turn on interpolate button
    plotVelocity(handles); %plot flow curves
    
    set(handles.loadCLpush,'Enable','on');
    set(handles.load2DPCbutton,'Enable','on');
    set(handles.pcPlanePopup,'Enable','on'); %make it so we can't select another plane
    set(handles.pcDatasetPopup,'Enable','on'); %make it so we can't select another plane
    set(handles.drawROIbutton,'Enable','on');
    set(handles.interpolatePopup,'Enable','on'); %set all possible interpolation types
    set(handles.errorBarRadio,'Enable','on');
    set(handles.pgShiftRadio,'Enable','on');

    guidata(hObject,handles);
    updatePCImages(handles); %update images (to remove green ROI circle)
    axes(handles.pcPlanePlot); %make sure we're still on PC plot
   

    % --- PLANE SLIDER - CALLBACK
function pcSlider_Callback(hObject, eventdata, handles)
    updatePCImages(handles); %update images if slider is moved

% --- PLANE SLIDER - CREATE FUNCTION
function pcSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
    
% --- MINIMUM CONTRAST VALUE BOX - CALLBACK  
    function minContrastBox_Callback(hObject, eventdata, handles)
    updatePCImages(handles)

% --- MINIMUM CONTRAST VALUE BOX - CREATE FUNCTION   
function minContrastBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- MAXIMUM CONTRAST VALUE BOX - CALLBACK  
function maxContrastBox_Callback(hObject, eventdata, handles)
    updatePCImages(handles)
    
% --- MAXIMUM CONTRAST VALUE BOX - CREATE FUNCTION
function maxContrastBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
    
%%%%%%%%%%%% VELOCITY PLOT %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- VELOCITY PLOT - CREATE FUNCTION
function velocityPlot_CreateFcn(hObject, eventdata, handles)
    

% --- INTERPOLATE POPUP - CALLBACK
function interpolatePopup_Callback(hObject, eventdata, handles)
    interp = get(handles.interpolatePopup,'Value');
    switch interp
        case 1
            handles.global.interpType = 'Gaussian'; %set global flag
        case 2
            handles.global.interpType = 'None';
    end 
 
    guidata(hObject, handles);
    plotVelocity(handles) %replot our velocity with interpolated data
    
    if handles.global.startAnalyzing %if we're already analyzing PWVs
        completeLoadingROI_Callback(hObject, eventdata, handles); %recompute PWVs with interpolated data
    end 

% --- INTERPOLATE POPUP - CREATE FUNCTION
function interpolatePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- ERROR BAR RADIO - CALLBACK
function errorBarRadio_Callback(hObject, eventdata, handles)
    plotVelocity(handles) %replot flow curves with errors bars
    
    
% --- PG SHIFT RADIO - CALLBACK
function pgShiftRadio_Callback(hObject, eventdata, handles)
    if get(handles.pgShiftRadio,'Value')
        handles.global.pgShift = 1;
    else
        handles.global.pgShift = 0;
    end 
    
    plotVelocity(handles) %replot flow curves with half cycle shift
    if handles.global.startAnalyzing %if we're already analyzing PWVs
        completeLoadingROI_Callback(hObject, eventdata, handles); %recompute PWVs with shifted curves
    end 
    guidata(hObject, handles);

   
    
% --- Executes on button press in completeLoadingROI.
function completeLoadingROI_Callback(hObject, eventdata, handles)
    handles.global.startAnalyzing = 1; %turn on flag to state that we are ready for PWV analysis
    handles.flow = organizeFlowInfo(handles);
    
    % COMPUTE TIME SHIFTS
    if ~isfield(handles.flow,'TTUpstroke')
        flow = computeTTs(handles.flow,handles.global);
        handles.flow = flow;
    end 
    
    % COMPUTE PWVs
    distance = [0 cumsum(handles.centerline.PlaneDistances)];
    TTPoint = [handles.flow.TTPoint]; 
    TTFoot = [handles.flow.TTFoot]; 
    TTUpstroke = [handles.flow.TTUpstroke]; 
    Xcorr = [handles.flow.Xcorr];

    numMethods = get(handles.ttpointRadio,'Value') + ...
        get(handles.ttfRadio,'Value') + ...
        get(handles.ttuRadio,'Value') + ...
        get(handles.xcorrRadio,'Value'); %add all PWV buttons turned on
    numCompares = numel(distance); %get number of time shift methods  
    average = zeros(1,numCompares); %initialize average timeshift array
    
    cla(handles.TimeVsDistance,'reset'); %reset PWV plot
    axes(handles.TimeVsDistance); hold on; %force axes to PWV plot
    xlabel('Distance (mm)'); ylabel('Time Shift (ms)'); %label axes
    sz = 30; %create marker sizes (filled in circles) of 30 pixels
    legendSet = {}; %initialize legend cell array
    if get(handles.ttpointRadio,'Value')
        scatter(distance,TTPoint,sz,'filled','MarkerFaceColor',[0.850 0.325 0.098]); %orange 
        for i=1:numCompares
            average(i) = average(i)+TTPoint(i); %add all distances for each ROI location
        end 
        legendSet{end+1} = 'TTPoint';
    end 
    if get(handles.ttfRadio,'Value')
        scatter(distance,TTFoot,sz,'filled','MarkerFaceColor',[0.494 0.184 0.556]); %purple 
        for i=1:numCompares
            average(i) = average(i)+TTFoot(i); %keep adding distances
        end 
        legendSet{end+1} = 'TTFoot';
    end 
    if get(handles.ttuRadio,'Value')
        scatter(distance,TTUpstroke,sz,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880]); %green
        for i=1:numCompares
            average(i) = average(i)+TTUpstroke(i);
        end 
        legendSet{end+1} = 'TTUpstroke';
    end 
    if get(handles.xcorrRadio,'Value')
        scatter(distance,Xcorr,sz,'filled','MarkerFaceColor',[0.6350 0.0780 0.1840]); %red 
        for i=1:numCompares 
            average(i) = average(i)+Xcorr(i);
        end 
        legendSet{end+1} = 'XCorr';
    end 

    D = max(distance);
    Tmin = min([TTPoint TTFoot TTUpstroke Xcorr]);
    Tmax = max([TTPoint TTFoot TTUpstroke Xcorr]);
    xlim([0 D+D*0.2]); ylim([Tmin-Tmin*0.02 Tmax+Tmax*0.02]);
        
    average = average./numMethods; %get average TT for each ROI
    scatter(distance,average,40,'black'); %open black circles (size 40 pixels)
    legendSet{end+1} = 'AVERAGE'; %add average to legend
    
    d = 0:round(D+D*0.2);
    % Here, we perform linear regression to find best fit line for each
    % time-to (TT) method. Note: lineFit(1)=slope; lineFit(2)=intercept

    %TTPoint Method
    PointFit = polyfit(distance,TTPoint,1); %fit line to TTPoint points
    linePoint = PointFit(1)*d + PointFit(2); %calculate line

    %TTFoot Method
    FootFit = polyfit(distance,TTFoot,1);
    lineFoot = FootFit(1)*d + FootFit(2);

    %TTUpstroke Method
    UpstrokeFit = polyfit(distance,TTUpstroke,1);
    lineUpstroke = UpstrokeFit(1)*d + UpstrokeFit(2);

    %Xcorr Method
    XcorrFit = polyfit(distance,Xcorr,1);
    lineXcorr = XcorrFit(1)*d + XcorrFit(2);

    %Average TT
    AverageFit = polyfit(distance,average,1);
    lineAverage= AverageFit(1)*d + AverageFit(2);

    PWVpoint = 1/PointFit(1); %PWV = 1/slope (mm/ms = m/s)
    PWVfoot = 1/FootFit(1);
    PWVupstroke = 1/UpstrokeFit(1);
    PWVxcorr = 1/XcorrFit(1);
    PWVaverage = 1/AverageFit(1);

    hold on; 
    if get(handles.ttpointRadio,'Value') %if our ttpoint button is on, plot average ttp
        plot(d,linePoint,':','LineWidth',0.2,'MarkerFaceColor',[0.8500 0.3250 0.0980]); %orange 
        set(handles.ttpointData,'String',[num2str(round(PWVpoint,2)) ' m/s']); %write out PWV value in text field
    end 
    if get(handles.ttfRadio,'Value')
        plot(d,lineFoot,':','LineWidth',0.2,'MarkerFaceColor',[0.4940 0.1840 0.5560]); %purple
        set(handles.ttfData,'String',[num2str(round(PWVfoot,2)) ' m/s']);
    end 
    if get(handles.ttuRadio,'Value')
        plot(d,lineUpstroke,':','LineWidth',0.2,'MarkerFaceColor',[0.4660 0.6740 0.1880]); %green
        set(handles.ttuData,'String',[num2str(round(PWVupstroke,2)) ' m/s']);
    end 
    if get(handles.xcorrRadio,'Value')
        plot(d,lineXcorr,':','LineWidth',0.2,'MarkerFaceColor',[0.6350 0.0780 0.1840]); %red
        set(handles.xcorrData,'String',[num2str(round(PWVxcorr,2)) ' m/s']);
    end 
    plot(d,lineAverage,'-k','LineWidth',0.2);
    legend(legendSet,'Location','northwest');
    hold off; 
    
    if numMethods>0 %if we have at least one ttbutton on
        set(handles.averageData,'String',[num2str(round(PWVaverage,2)) ' m/s']); %set PWV text field
    else %if have not methods selected (all ttbuttons are off)
        set(handles.averageData,'String','0 m/s'); %set PWV text field to 0 m/s
    end
    
    set(handles.ttpointRadio,'Enable','on');
    set(handles.ttuRadio,'Enable','on');
    set(handles.ttfRadio,'Enable','on');
    set(handles.xcorrRadio,'Enable','on');
    set(handles.exportAnalysisButton,'Enable','on');
    guidata(hObject, handles);
    
    

%%%%%%%%%%%% PWV PLOT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% --- PWV PLOT - CREATE FUNCTION
function TimeVsDistance_CreateFcn(hObject, eventdata, handles)


% --- EXPORT ANALYSIS - CALLBACK
function exportAnalysisButton_Callback(hObject, eventdata, handles)
    interpTypes = {'None','Gaussian'};
    for t=1:length(interpTypes) %export data for each interp type
        handles.global.interpType = interpTypes{t}; %switch interp type
        flow = computeTTs(handles.flow,handles.global); %recalculate flow struct for each interp
        numROIs = numel(flow);
        
        distance = [0 cumsum(handles.centerline.PlaneDistances)];
        ttpoint = [handles.flow.TTPoint]; 
        ttfoot = [handles.flow.TTFoot];
        ttupstroke = [handles.flow.TTUpstroke]; 
        xcorr = [handles.flow.Xcorr];
        ttaverage = mean([ttfoot; ttpoint; ttupstroke; xcorr],1); %get average time shift
        
        Distances = abs(diff([distance 0])); %add zero at end to get full distance
        TTPoint = [diff(ttpoint) sum(diff(ttpoint))];
        TTFoot = [diff(ttfoot) sum(diff(ttfoot))];
        TTUpstroke = [diff(ttupstroke) sum(diff(ttupstroke))];
        Xcorr = [diff(xcorr) sum(diff(xcorr))];
        TTaverage = [diff(ttaverage) sum(diff(ttaverage))];
        
        % Make first column for excel file
        numCompares = length(Distances); %get number of PWV measurements
        for d=1:(numCompares-1)
            PLANES{d} = [flow(d).Name ' --> ' flow(d+1).Name]; %get names for Excel
        end 
        PLANES{numCompares} = [flow(1).Name ' --> ' flow(end).Name]; %get names for Excel
        PLANES{numCompares+1} = 'FIT_PWV'; %These rows be for PWV fit parameters
        PLANES{numCompares+2} = 'm (slope)';
        PLANES{numCompares+3} = 'b (y-intercept)';
        PLANES{numCompares+4} = 'R^2';

        PWV_Point = NaN(numCompares+4,1); %add dummy rows to match PLANES size
        PWV_Foot = NaN(numCompares+4,1);
        PWV_Upstroke = NaN(numCompares+4,1);
        PWV_Xcorr = NaN(numCompares+4,1);
        PWV_Average = zeros(1,numCompares);
        
        for i=1:numCompares
            PWV_Point(i) = Distances(i)/TTPoint(i); %calculate pointwise PWVs (not fits)
            PWV_Foot(i) = Distances(i)/TTFoot(i);
            PWV_Upstroke(i) = Distances(i)/TTUpstroke(i);
            PWV_Xcorr(i) = Distances(i)/Xcorr(i);
            PWV_Average(i) = (PWV_Point(i)+PWV_Foot(i)+PWV_Upstroke(i)+PWV_Xcorr(i))/4; %get simple average of all PWVs
        end 
        
        [linePointFit,S] = polyfit(distance,ttpoint,1); %get linear regression fit for all points 
        PWV_Point(numCompares+1) = 1/linePointFit(1); %calculate PWV (=1/slope)
        PWV_Point(numCompares+2) = linePointFit(1); %get slope
        PWV_Point(numCompares+3) = linePointFit(2); %get y-intercept
        PWV_Point(numCompares+4) = (1 - (S.normr/norm(ttpoint - mean(ttpoint)))^2); %R^2

        [lineFootFit,S] = polyfit(distance,ttfoot,1);
        PWV_Foot(numCompares+1) = 1/lineFootFit(1);
        PWV_Foot(numCompares+2) = lineFootFit(1);
        PWV_Foot(numCompares+3) = lineFootFit(2);
        PWV_Foot(numCompares+4) = (1 - (S.normr/norm(ttfoot - mean(ttfoot)))^2);

        [lineUpstrokeFit,~] = polyfit(distance,ttupstroke,1);
        PWV_Upstroke(numCompares+1) = 1/lineUpstrokeFit(1);
        PWV_Upstroke(numCompares+2) = lineUpstrokeFit(1);
        PWV_Upstroke(numCompares+3) = lineUpstrokeFit(2);
        PWV_Upstroke(numCompares+4) = (1 - (S.normr/norm(ttupstroke - mean(ttupstroke)))^2);

        [lineXcorrFit,~] = polyfit(distance,xcorr,1);
        PWV_Xcorr(numCompares+1) = 1/lineXcorrFit(1);
        PWV_Xcorr(numCompares+2) = lineXcorrFit(1);
        PWV_Xcorr(numCompares+3) = lineXcorrFit(2);
        PWV_Xcorr(numCompares+4) = (1 - (S.normr/norm(xcorr - mean(xcorr)))^2);

        [lineAverageFit,~] = polyfit(distance,ttaverage,1);
        PWV_Average(numCompares+1) = 1/lineAverageFit(1);
        PWV_Average(numCompares+2) = lineAverageFit(1);
        PWV_Average(numCompares+3) = lineAverageFit(2);
        PWV_Average(numCompares+4) = (1 - (S.normr/norm(ttaverage - mean(ttaverage)))^2);
        PWV_Average = PWV_Average';
        
        Distances(numCompares+1:numCompares+4) = NaN;
        TTPoint(numCompares+1:numCompares+4) = NaN;
        TTFoot(numCompares+1:numCompares+4) = NaN;
        TTUpstroke(numCompares+1:numCompares+4) = NaN;
        Xcorr(numCompares+1:numCompares+4) = NaN;
        TTaverage(numCompares+1:numCompares+4) = NaN;
        
        PLANES = PLANES'; %Needed for excel, won't save name if ' in table call
        Distances = Distances';
        TTPoint = TTPoint';
        TTFoot = TTFoot';
        TTUpstroke = TTUpstroke';
        Xcorr = Xcorr';
        TTaverage = TTaverage';
        
        % Make table for writing excel file
        pwvTable = table(PLANES,Distances,TTPoint,TTFoot,TTUpstroke,Xcorr,TTaverage,PWV_Point,PWV_Foot,PWV_Upstroke,PWV_Xcorr,PWV_Average);
        baseDir = handles.global.homeDir; %rejoin string to get name of folder one up from plane data
        date = datestr(now); %get current date/time
        chopDate = [date(1:2) '-' date(4:6) '-' date(10:11) '-' date(13:14) date(16:17)]; %chop date up
        if ~exist([baseDir filesep 'DataAnalysis'],'dir') %if directory doesn't exist
            mkdir([baseDir filesep 'DataAnalysis']); %make it
            cd([baseDir filesep 'DataAnalysis']); %go to it
            writetable(pwvTable,['Summary_' chopDate '.xlsx'],'FileType','spreadsheet','Sheet',['Interpolation - ' interpTypes{t}]); %write excel sheet for each interp
            flow = handles.flow; %make variables for saving
            if strcmp(handles.global.interpType,'Gaussian')
                saveTTplots(handles,flow);
            end 
            save('flow.mat','flow')
            save('pwvTable.mat','pwvTable')
        else %or if the directory already exists
            cd([baseDir filesep 'DataAnalysis']); %go to it
            writetable(pwvTable,['Summary_' chopDate '.xlsx'],'FileType','spreadsheet','Sheet',['Interpolation - ' interpTypes{t}]);
            flow = handles.flow;
            if strcmp(handles.global.interpType,'Gaussian')
                saveTTplots(handles,flow);
            end 
            save('flow.mat','flow')
            save('pwvTable.mat','pwvTable')
        end 
        cd(handles.global.homeDir); %go back home  
        clear PLANES %need to do this because PLANES will keep getting transposed
    end 
    
    cd([baseDir filesep 'DataAnalysis']); %go to it
    frame = getframe(handles.TimeVsDistance); %get snapshot of PWV plot 
    imwrite(frame2im(frame),'PWVanalysisPlot.png'); %write out to PNG
    pcDatasets = handles.pcDatasets;
    save('pcDatasets.mat','pcDatasets');
    cd(handles.global.homeDir); %go back home  
    set(handles.exportDone,'String','Export Completed!');

    guidata(hObject, handles);
    
   

% --- TTPoint READOUT - CREATE FUNCTION
function ttpointData_CreateFcn(hObject, eventdata, handles)
% --- TTPoint RADIO - CALLBACK
function ttpointRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttpointRadio,'Value') %if we're turned off
        set(handles.ttpointData,'String',' '); %don't display PWV
    end 
    if handles.global.startAnalyzing %if we're analyzing PWVs
        completeLoadingROI_Callback(hObject, eventdata, handles); %reanalyze without TTpoint 
    end 

% --- TTUpstroke READOUT - CREATE FUNCTION
function ttuData_CreateFcn(hObject, eventdata, handles)
% --- TTUpstroke RADIO - CALLBACK
function ttuRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttuRadio,'Value')
        set(handles.ttuData,'String',' ');
    end 
    if handles.global.startAnalyzing
        completeLoadingROI_Callback(hObject, eventdata, handles);
    end 

% --- TTFoot READOUT - CREATE FUNCTION
function ttfData_CreateFcn(hObject, eventdata, handles)
% --- TTFoot RADIO - CALLBACK
function ttfRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttfRadio,'Value')
        set(handles.ttfData,'String',' ');
    end 
    if handles.global.startAnalyzing
        completeLoadingROI_Callback(hObject, eventdata, handles);
    end 

% --- Xcorr READOUT - CREATE FUNCTION
function xcorrData_CreateFcn(hObject, eventdata, handles)
% --- Xcorr RADIO - CALLBACK
function xcorrRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.xcorrRadio,'Value')
        set(handles.xcorrData,'String',' ');
    end 
    if handles.global.startAnalyzing
        completeLoadingROI_Callback(hObject, eventdata, handles);
    end 

% --- AVERAGE READOUT - CREATE FUNCTION
function averageData_CreateFcn(hObject, eventdata, handles)





%%%%%%%%%%%% MY FUNCTIONS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% --- Update Images in PLANE PLOT
function updatePCImages(handles)
    axes(handles.pcPlanePlot); %force axes to PC plot
    
    planeNum = get(handles.pcPlanePopup,'Value'); % get current plane index (eg AAo)
    datasetNum = get(handles.pcDatasetPopup,'Value'); %get current dataset (eg MAG)
    
    dataset = handles.pcDatasets;
    if isstruct(dataset(planeNum).Images) %if Data is a structure
        imageSet = struct2cell(dataset(planeNum).Images); %convert from struct to cell
    else  %or if we already have an array
        imageSet = dataset(planeNum).Images; %just change the name
    end 

    if iscell(imageSet) %if our set of images are contained in a cell
        images = imageSet(datasetNum); %pull images for current dataset
        images = cell2mat(images); %turn to matrix
    else 
        images = imageSet; %otherwise, do nothing
    end 
    
    if ndims(images)<3 %if we are dealing with time-averaged images (ndim=2)
        maxSize = max(size(images,1),size(images,2));
        steps = [1 maxSize]; %set our slider as wide as possible so we can't slide
        set(handles.pcSlider,'SliderStep', steps);
        slice = images; %change name for consistency (see below)
    else %if we are dealing with time-resolved images (ndim=3)
        dim3size = size(images,3); %get the size of the third dimension
        steps = [1/(dim3size-1) 10/(dim3size-1)]; %set so one 'slide' moves to the next slice exactly
        set(handles.pcSlider,'SliderStep', steps);
        sliceNum = 1+round( get(handles.pcSlider,'Value').*(dim3size-1) ); %get slice number from slider
        slice = images(:,:,sliceNum); %pull slice from images
    end 
    
    minc = str2double(get(handles.minContrastBox,'String'));
    maxc = str2double(get(handles.maxContrastBox,'String'));

    if ~isempty(handles.pcDatasets(planeNum).ROI) %if we have an ROI placed
        hold on %make it so we can slide while keeping the ROI on the figure
        imshow(rescale(slice),[minc maxc]);
    else 
        cla(handles.pcPlanePlot,'reset') %otherwise, reset the plot
        imshow(rescale(slice),[minc maxc]) %then show the image
    end 
    
    
    
% --- "Time to" calculations (TTPoint, TTUpstroke, TTFoot, Xcorr)
function flow = computeTTs(flow,globals)
    pgShift = globals.pgShift;
    numROIs = globals.totalROIs;
    for i=1:numROIs %for each ROI
        if globals.startAnalyzing %if we're analyzing (put here because we call computePWV in interpolatePopup)
            switch globals.interpType  %find the appropriate interp data
                case 'None'
                    flowTemp = flow(i).Uninterp.meanROI;
                case 'Gaussian'
                    flowTemp = flow(i).Gaussian.meanROI;  
            end 
        else %else, we just use uninterpolated data
            flowTemp = flow(i).Uninterp.meanROI;
        end 
        
        if mean(flowTemp)<0 %if our flow curve is mainly negative
            flowTemp = -1*flowTemp; %flip it upside down so we can find time shifts
        end 
        
        if pgShift %if we want to shift our waveform over
            flowTemp = circshift(flowTemp,round(length(flowTemp)/3)); %shift by half cycle
        end 
                
%         times = flow(i).Uninterp.times; %get time frames (ms)
        times = flow(1).Uninterp.times;
        timeres = times(2)-times(1); %temporal resolution (ms)
%         curvePoints(i).times = times;
%         curvePoints(i).timeres = timeres;
        
        [maxPeakVel,maxPeakVelIdx] = max(flowTemp); %find max velocity value and its location
        upstroke = flowTemp(1:maxPeakVelIdx); %define 'upstroke' region of the flow curve
        curvePoints(i).maxPeakVelIdx = maxPeakVelIdx; %add max point to curvePoints struct
        curvePoints(i).maxPeakVel = maxPeakVel; %add max velocity to curvePoints struct
        
        [~,EightyPointIdx] = min(abs(upstroke-0.8*maxPeakVel)); %get point at 70% max peak
        curvePoints(i).EightyPointIdx = EightyPointIdx; %add to curvePoints struct
        curvePoints(i).EightyPoint = flowTemp(EightyPointIdx); %add 70% flow value to curvePoints
        
        [~,TwentyPointIdx] = min(abs(upstroke-0.2*maxPeakVel)); %get point at 30% max peak
        curvePoints(i).TwentyPointIdx = TwentyPointIdx; %add to curvePoints struct
        curvePoints(i).TwentyPoint = flowTemp(TwentyPointIdx);
        
        [~,FiftyPointIdx] = min(abs(upstroke-0.5*maxPeakVel)); %get point at 50% max peak    
        curvePoints(i).FiftyPointIdx = FiftyPointIdx;
        curvePoints(i).FiftyPoint = flowTemp(FiftyPointIdx);

        flows(i,:) = normalize(flowTemp,'range'); %normalize curve from here on
    end      
    
    % TTPoint - time to point calculation
    for i=1:numROIs
        tp = times(curvePoints(i).FiftyPointIdx); %get time at 50% peak (ms) 
        flow(i).TTPoint = tp; %add to flow struct
    end    
    
    % TTUpstroke - time to upstroke calculation
    for i=1:numROIs
        [sigmoid,sigTimes,tp] = sigFit(flows(i,:),times); %see sigFit function below
        flow(i).TTUpstroke = tp; %add to flow struct
        flow(i).SigmoidFit = sigmoid;
        flow(i).SigmoidTimes = sigTimes;
    end  
  
    % TTFoot - time to foot calculation
    for i=1:numROIs
        leftIdx = curvePoints(i).TwentyPointIdx;
        rightIdx = curvePoints(i).EightyPointIdx;
        flowSeg = flows(i,leftIdx:rightIdx);
        timeSeg = times(leftIdx:rightIdx);
        p1 = polyfit(timeSeg,flowSeg,1);
        tp = -(p1(2)/p1(1));
        flow(i).TTFoot = tp; %add to flow struct
        flow(i).P1 = p1;
    end  
  
    % XCorr - cross correlation calculation
    %since we need 2 curves to get shift, make 1st point = to 1st TTF
    flow(1).Xcorr = flow(1).TTFoot;
    for i=2:numROIs
        flow1 = flows(i-1,:);
        flow2 = flows(i,:);
        [Xcorrs,lags] = xcorr(flow2,flow1,'normalized'); %perform cross correlation between flow curves
        [~,maxXcorrIdx] = max(Xcorrs); %get index of max Xcorr value
        shift = lags(maxXcorrIdx); %find time lag of Xcorr peak
        flow(i).Xcorr = flow(i-1).Xcorr + timeres*shift; %cumulative sum from last time
    end  
    
    
 
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
    

    
% --- Plot Velocities    
function plotVelocity(handles)
    cla(handles.velocityPlot,'reset'); %reset axes
    axes(handles.velocityPlot); %make sure we plot on the right axis
    
    times = handles.pcDatasets(1).ROIdataUninterp(1).times;
    plot(times, zeros(1,length(times)) ,'Color','black','LineWidth',1.5); %line of y=0 (for visual reference)
    % Note that above we assume same time scale for each plane (MR scan)
    xlim([min(times) max(times)]);
    legendSet = {'Baseline'}; %add baseline to legend names
    for i=1:numel(handles.pcDatasets) %for all planes
        if isstruct(handles.pcDatasets(i).ROIdataUninterp) %if we've made ROI data for this dataset
            for j=1:length(handles.pcDatasets(i).ROIdataUninterp) %for each ROI
                legendSet{end+1} = handles.pcDatasets(i).ROIinfo(j).Name; %add name of ROI to list
                switch handles.global.interpType %check what interpolation we are using
                    case 'None'
                        velocity = handles.pcDatasets(i).ROIdataUninterp(j).meanROI; %grab mean velocity
                        %times = handles.pcDatasets(i).ROIdataUninterp(j).times;
                        stdv = handles.pcDatasets(i).ROIdataUninterp(j).stdvROI; %grab stdv of velocity
                    case 'Gaussian'
                        velocity = handles.pcDatasets(i).ROIdataGaussian(j).meanROI;
                        %times = handles.pcDatasets(i).ROIdataUninterp(j).times;
                        stdv = handles.pcDatasets(i).ROIdataGaussian(j).stdvROI;
                end 
                
                if get(handles.pgShiftRadio,'Value') %check if we need to shift waveform
                    velocity = circshift(velocity,round(length(velocity)/3)); %shift by half cycle
                    stdv = circshift(stdv,round(length(stdv)/3));
                end 

                if mean(velocity)<0 %if we are mainly negative velocities (as in descending aorta)
                    if get(handles.errorBarRadio,'Value') %and if we want to show error bars
                        hold on; errorbar(times,-1*velocity,stdv); %plot inverted velocity and errors
                    else
                        hold on; plot(times,-1*velocity); %else, plot inverted velocity
                    end 
                else %otherwise, don't invert velocity (as in ascending aorta)
                    if get(handles.errorBarRadio,'Value')
                        hold on; errorbar(times,velocity,stdv);
                    else
                        hold on; plot(times,velocity);
                    end
                end 
            end 
        end 
    end 
    legend(legendSet); hold off
    xlabel('Time (ms)'); ylabel('Mean Velocity in ROI (mm/s)'); %set axes labels
    %xlim([times(1),times(end)]); %chop limits to make curve full width
    

    
% --- Sigmoid Fit Function
function [sigmoid,t,t1] = sigFit(meanROI,times)
%%% See the following article by Anas Dogui in JMRI:
% Measurement of Aortic Arch Pulse Wave Velocity in Cardiovascular MR:
% Comparison of Transit Time Estimators and Description of a New Approach

    [~,peak] = max(meanROI); %find max
    upslope = meanROI(1:peak); %find upslope region of flow curve
    t0 = times(1:peak);
    t = linspace(1,times(peak),1000); %interpolate even more
    upslope = interp1(t0,upslope,t); %interpolate upslope
    upslope = rescale(upslope); %normalize from 0 to 1
    [~,MIN] = min(upslope);
    upslope = upslope(MIN:end);
    t = t(MIN:end);
    dt = t(2)-t(1); %new temporal resolution (=0.1)
    midpoint = round(length(upslope)/2);
    
    % c1 = b, c2 = a, c3 = x0, c4 = dx
    % Note that we could assume the equation e^t/(1+e^(t-t0)) since c1=1 and
    % c2=0. However, will keep the same as the Dogui paper.
    sigmoidModel = @(c) c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((t-c(3))./c(4)) ) - upslope;
    c0 = [0,1,t(midpoint),dt/2]; %initial params for upslope region
    opts = optimset('Display', 'off'); %turn off display output
    c = lsqnonlin(sigmoidModel,c0,[],[],opts); %get nonlinear LSQ solution
    sigmoid = c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((t-c(3))./c(4)) ); %calculate our sigmoid fit with our new params

    dy = diff(sigmoid,1);
    dy(end) = [];
    dx = diff(t,1);
    dx(end) = [];
    ddy = diff(sigmoid,2);
    curvature = ddy.*dx./(dx.^2 + dy.^2).^(3/2);
    [~,tIdx] = max(curvature);
    t1 = t(tIdx);

    
    
% --- Condense and organize flow data obtained from ROIs
function flow = organizeFlowInfo(handles)
% This function is designed to pull apart handles.pcDatasets. It is
% difficult to do analysis on the handles structure because some slices
% have 2 ROIs. It is much easier to have a structure that pulls out each
% ROI. It only adds a bit of memory since we aren't saving the raw images.
    count = 1; %overall iterator
    for i=1:numel(handles.pcDatasets) %for each PC dataset
        if isstruct(handles.pcDatasets(i).ROIdataUninterp) %do we have data?
            if length(handles.pcDatasets(i).ROIdataUninterp)==1 %if we just have one ROI
                flow(count).Name = handles.pcDatasets(i).ROIinfo.Name; %pull only relevant info for PWV calcs
                flow(count).ROIinfo = handles.pcDatasets(i).ROIinfo;
                flow(count).HeaderInfo = handles.pcDatasets(i).Info;
                flow(count).Uninterp = handles.pcDatasets(i).ROIdataUninterp;
                flow(count).Gaussian = handles.pcDatasets(i).ROIdataGaussian;
                flow(count).pcDatasetREF = [i,1];
                count = count+1;
            else 
                for j=1:length(handles.pcDatasets(i).ROIdataUninterp) %if we have more than one ROI
                    flow(count).Name = handles.pcDatasets(i).ROIinfo(j).Name; %parse into individual ROIs in flow struct
                    flow(count).ROIinfo = handles.pcDatasets(i).ROIinfo(j);
                    flow(count).HeaderInfo = handles.pcDatasets(i).Info;
                    flow(count).Uninterp = handles.pcDatasets(i).ROIdataUninterp(j);
                    flow(count).Gaussian = handles.pcDatasets(i).ROIdataGaussian(j);
                    flow(count).pcDatasetREF = [i,j];
                    count = count+1;
                end 
            end 
        end 
    end 

 function saveTTplots(handles,flow)   
    times = handles.pcDatasets(1).ROIdataUninterp(1).times;
    figure('units','normalized','outerposition',[0 0 1 1]); 
    plot(times, zeros(1,length(times)) ,'Color','black','LineWidth',1.5); %line of y=0 (for visual reference)
    % Note that above we assume same time scale for each plane (MR scan)
    xlim([min(times) max(times)]);
    legendSet = {'Baseline'}; %add baseline to legend names
    count = 1;
    for i=1:numel(handles.pcDatasets) %for all planes
        if isstruct(handles.pcDatasets(i).ROIdataUninterp) %if we've made ROI data for this dataset
            for j=1:length(handles.pcDatasets(i).ROIdataUninterp) %for each ROI
                legendSet{end+1} = handles.pcDatasets(i).ROIinfo(j).Name; %add name of ROI to list
                switch handles.global.interpType %check what interpolation we are using
                    case 'None'
                        vTemp = handles.pcDatasets(i).ROIdataUninterp(j).meanROI; %grab mean velocity
                    case 'Gaussian'
                        vTemp = handles.pcDatasets(i).ROIdataGaussian(j).meanROI;
                end 
                
                if get(handles.pgShiftRadio,'Value') %check if we need to shift waveform
                    vTemp = circshift(vTemp,round(length(vTemp)/3)); %shift by half cycle
                end 

                if mean(vTemp)<0
                    vTemp = -1*vTemp;
                end
                
                hold on; plot(times,vTemp);
                velocity(count,:) = vTemp;
                count = count+1;
            end 
        end 
    end 
    legend(legendSet,'Location','northeastoutside','AutoUpdate','off');
    xlabel('Time (ms)'); 
    mx = max(velocity(:));
    mn = min(velocity(:));
    ylim([(mn - mn*0.05) (mx + mx*0.05)]);
    ylabel('Mean Velocity in ROI (mm/s)'); %set axes labels
    ttp = [flow.TTPoint];
    for i=1:length(ttp)
        [~,idx] = min( abs(times-ttp(i)) );
        h(i) = scatter(times(idx),velocity(i,idx),'k','filled');
        h2(i) = scatter(times(idx),0,'k');
        h3(i) = plot([times(idx) times(idx)],[velocity(i,idx),0],':k');
    end 
    title('Time-To-Point');
    frame = getframe(gcf);
    imwrite(frame2im(frame),[pwd filesep 'TTPoint.png']);
    
    set(h,'Visible','off'); set(h2,'Visible','off'); set(h3,'Visible','off')
    
    ttf = [flow.TTFoot];
    for i=1:length(ttf)
        [~,idx] = min( abs(times-ttf(i)) );
        [scale,idx2] = max(velocity(i,:));
        P1 = flow(i).P1;
        x = times(idx:idx2);
        y = scale*(P1(1)*x + P1(2));
        h(i) = plot(x,y,'k');
        h2(i) = scatter(times(idx),0,'k');
    end 
    title('Time-To-Foot');
    frame = getframe(gcf);
    imwrite(frame2im(frame),[pwd filesep 'TTFoot.png']);
    
    set(h,'Visible','off'); set(h2,'Visible','off');
        
    ttu = [flow.TTUpstroke];
    for i=1:length(ttu)
        [~,idx] = min( abs(times-ttu(i)) );
        [scale,~] = max(velocity(i,:));
        h(i) = scatter(times(idx),velocity(i,idx),'k','filled');
        h2(i) = scatter(times(idx),0,'k');
        h3(i) = plot([times(idx) times(idx)],[velocity(i,idx),0],':k');
        h4(i) = plot(flow(i).SigmoidTimes,scale*flow(i).SigmoidFit,'.k');
    end 
    title('Time-To-Upstroke');
    frame = getframe(gcf);
    imwrite(frame2im(frame),[pwd filesep 'TTUpstroke.png']);
    
    set(h,'Visible','off'); set(h2,'Visible','off'); set(h3,'Visible','off'); set(h4,'Visible','off');
    
    xcorr = [flow.Xcorr];
    for i=1:length(xcorr)
        [~,idx] = min( abs(times-xcorr(i)) );
        h(i) = scatter(times(idx),velocity(i,idx),'k','filled');
        h2(i) = scatter(times(idx),0,'k');
        h3(i) = plot([times(idx) times(idx)],[velocity(i,idx),0],':k');
    end 
    title('Cross Correlation Time Lag');
    frame = getframe(gcf);
    imwrite(frame2im(frame),[pwd filesep 'Xcorr.png']);
    
    set(h,'Visible','off'); set(h2,'Visible','off'); set(h3,'Visible','off')
    
    hold off; close(gcf);   
    
    
    
% Load Dat files
function v = load_dat(name, res)
    [fid,errmsg]= fopen(name,'r');
    if fid < 0  %if name does not exist in directory
        set(handles.MessageBar,'String',['Error Opening Data : ',errmsg]);
    end

    % Reads in as short, reshapes by image res.
    v = reshape(fread(fid,'short=>single'),res);
    fclose(fid);

basedir = pwd;
now = datestr(now,'mm/dd/yyyy');

%% Cartesian
if exist('PWV_2DPC_Analysis\DataAnalysis_Cartesian','dir')
    load('PWV_2DPC_Analysis\DataAnalysis_Cartesian\flow.mat');
    roi_cart = length(flow);
    info = flow(1).HeaderInfo;
    dt = datetime([info.AcquisitionDate ' ' info.AcquisitionTime],'InputFormat','yyyyMMdd HHmmss');
    scan_date = datestr(dt);
    mri_loc = info.InstitutionName;
    rrInt = info.NominalInterval;
    frames = info.CardiacNumberOfImages;
    TR_AAo_CART = rrInt/frames;
    TR_AbdAo_CART = NaN;
    plane_dx = NaN;

    load('PWV_2DPC_Analysis\DataAnalysis_Cartesian\pwvTable.mat');
    dx1_2 = pwvTable{1,2};
    dx2_3 = NaN;
    dx1_3 = NaN;
    ttp_pwv_CART = pwvTable{3,8};
    ttf_pwv_CART = pwvTable{3,9};
    ttu_pwv_CART = pwvTable{3,10};
    xcor_pwv_CART = pwvTable{3,11};
else
    scan_date = NaN;
    mri_loc = NaN;
    TR_AAo_CART = NaN;
    TR_AbdAo_CART = NaN;
    plane_dx = NaN;
    dx1_2 = NaN;
    dx2_3 = NaN;
    dx1_3 = NaN;
    ttp_pwv_CART = NaN;
    ttf_pwv_CART = NaN;
    ttu_pwv_CART = NaN;
    xcor_pwv_CART = NaN;
end 
% roi_rad = NaN;
% params = {roi_rad roi_cart now scan_date mri_loc dx1_2 dx2_3 dx1_3 plane_dx ...
%     ttp_pwv_CART ttf_pwv_CART ttu_pwv_CART xcor_pwv_CART TR_AAo_CART TR_AbdAo_CART};
% clear basedir dx* flow frames info mri_loc now plane_dx pwvTable roi* rrInt scan_date
% clear TR* tt* xcor* dt

%% Radial Low Resolution
if exist('PWV_2DPC_Analysis\DataAnalysis_Radial_LowRes','dir')
    load('PWV_2DPC_Analysis\DataAnalysis_Radial_LowRes\flow.mat');
    roi_rad = length(flow);
    info = flow(1).HeaderInfo;
    TR_AAo_radLR = info.timeres;
    TR_AbdAo_radLR = NaN;

    load('PWV_2DPC_Analysis\DataAnalysis_Radial_LowRes\pwvTable.mat');
    ttp_pwv_radLR = pwvTable{3,8};
    ttf_pwv_radLR = pwvTable{3,9};
    ttu_pwv_radLR = pwvTable{3,10};
    xcor_pwv_radLR = pwvTable{3,11};
else
    TR_AAo_radLR = NaN;
    TR_AbdAo_radLR = NaN;
    ttp_pwv_radLR = NaN;
    ttf_pwv_radLR = NaN;
    ttu_pwv_radLR = NaN;
    xcor_pwv_radLR = NaN;
end 

%% Radial High Resolution (Local Low Rank Reconstruction)
if exist('PWV_2DPC_Analysis\DataAnalysis_Radial_HighRes','dir')
    load('PWV_2DPC_Analysis\DataAnalysis_Radial_HighRes\flow.mat');
    info = flow(1).HeaderInfo;
    TR_AAo_radHR = info.timeres;
    TR_AbdAo_radHR = NaN;

    load('PWV_2DPC_Analysis\DataAnalysis_Radial_HighRes\pwvTable.mat');
    ttp_pwv_radHR = pwvTable{3,8};
    ttf_pwv_radHR = pwvTable{3,9};
    ttu_pwv_radHR = pwvTable{3,10};
    xcor_pwv_radHR = pwvTable{3,11};
else
    TR_AAo_radHR = NaN;
    TR_AbdAo_radHR = NaN;
    ttp_pwv_radHR = NaN;
    ttf_pwv_radHR = NaN;
    ttu_pwv_radHR = NaN;
    xcor_pwv_radHR = NaN;
end 

params = {roi_rad roi_cart now scan_date mri_loc dx1_2 dx2_3 dx1_3 plane_dx ...
    ttp_pwv_CART ttf_pwv_CART ttu_pwv_CART xcor_pwv_CART TR_AAo_CART TR_AbdAo_CART ...
    ttp_pwv_radLR ttf_pwv_radLR ttu_pwv_radLR xcor_pwv_radLR TR_AAo_radLR TR_AbdAo_radLR ...
    ttp_pwv_radHR ttf_pwv_radHR ttu_pwv_radHR xcor_pwv_radHR TR_AAo_radHR TR_AbdAo_radHR};
clear basedir dx* flow frames info mri_loc now plane_dx pwvTable roi* rrInt scan_date
clear TR* tt* xcor* dt 
basedir = pwd;
cd('PWV_2DPC_Analysis')

%% Cartesian
if exist('DataAnalysis_Cartesian','dir')
    load('DataAnalysis_Cartesian\flow.mat');
    info = flow(1).HeaderInfo;
    scan_date = info.FileModDate;
    mri_loc = info.InstitutionName;
    rrInt = info.NominalInterval;
    frames = info.CardiacNumberOfImages;
    TR_AAo_CART = rrInt/frames;
    dx1 = info.ImagePositionPatient(3);
    
    info = flow(3).HeaderInfo;
    rrInt = info.NominalInterval;
    frames = info.CardiacNumberOfImages;
    TR_AbdAo_CART = rrInt/frames;
    dx2 = info.ImagePositionPatient(3);
    plane_dx = abs(dx1-dx2); 

    load('DataAnalysis_Cartesian\pwvTable.mat');
    dx1_2 = pwvTable{1,2};
    dx2_3 = pwvTable{2,2};
    dx1_3 = pwvTable{3,2};
    ttp_pwv_CART = pwvTable{4,8};
    ttf_pwv_CART = pwvTable{4,9};
    ttu_pwv_CART = pwvTable{4,10};
    xcor_pwv_CART = pwvTable{4,11};
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

%% Radial Low Resolution
if exist('DataAnalysis_Radial_LowRes','dir')
    load('DataAnalysis_Radial_LowRes\flow.mat');
    info = flow(1).HeaderInfo;
    TR_AAo_radLR = info.timeres;

    info = flow(3).HeaderInfo;
    TR_AbdAo_radLR = info.timeres;

    load('DataAnalysis_Radial_LowRes\pwvTable.mat');
    ttp_pwv_radLR = pwvTable{4,8};
    ttf_pwv_radLR = pwvTable{4,9};
    ttu_pwv_radLR = pwvTable{4,10};
    xcor_pwv_radLR = pwvTable{4,11};
else
    TR_AAo_radLR = NaN;
    TR_AbdAo_radLR = NaN;
    ttp_pwv_radLR = NaN;
    ttf_pwv_radLR = NaN;
    ttu_pwv_radLR = NaN;
    xcor_pwv_radLR = NaN;
end 

%% Radial Low Resolution - Expiration
if exist('DataAnalysis_Radial_Expiration','dir')
    load('DataAnalysis_Radial_Expiration\flow.mat');
    info = flow(1).HeaderInfo;
    TR_AAo_radEXP = info.timeres;

    info = flow(3).HeaderInfo;
    TR_AbdAo_radEXP = info.timeres;

    load('DataAnalysis_Radial_Expration\pwvTable.mat');
    ttp_pwv_radEXP = pwvTable{4,8};
    ttf_pwv_radEXP = pwvTable{4,9};
    ttu_pwv_radEXP = pwvTable{4,10};
    xcor_pwv_radEXP = pwvTable{4,11};
else
    TR_AAo_radEXP = NaN;
    TR_AbdAo_radEXP = NaN;
    ttp_pwv_radEXP = NaN;
    ttf_pwv_radEXP = NaN;
    ttu_pwv_radEXP = NaN;
    xcor_pwv_radEXP = NaN;
end 

%% Radial High Resolution (Local Low Rank Reconstruction)
if exist('DataAnalysis_Radial_HighRes','dir')
    load('DataAnalysis_Radial_HighRes\flow.mat');
    info = flow(1).HeaderInfo;
    TR_AAo_radHR = info.timeres;

    info = flow(3).HeaderInfo;
    TR_AbdAo_radHR = info.timeres;

    load('DataAnalysis_Radial_HighRes\pwvTable.mat');
    ttp_pwv_radHR = pwvTable{4,8};
    ttf_pwv_radHR = pwvTable{4,9};
    ttu_pwv_radHR = pwvTable{4,10};
    xcor_pwv_radHR = pwvTable{4,11};
else
    TR_AAo_radHR = NaN;
    TR_AbdAo_radHR = NaN;
    ttp_pwv_radHR = NaN;
    ttf_pwv_radHR = NaN;
    ttu_pwv_radHR = NaN;
    xcor_pwv_radHR = NaN;
end 

params = {scan_date mri_loc dx1_2 dx2_3 dx1_3 plane_dx ...
    ttp_pwv_CART ttf_pwv_CART ttu_pwv_CART xcor_pwv_CART TR_AAo_CART TR_AbdAo_CART ...
    ttp_pwv_radLR ttf_pwv_radLR ttu_pwv_radLR xcor_pwv_radLR TR_AAo_radLR TR_AbdAo_radLR ...
    ttp_pwv_radEXP ttf_pwv_radEXP ttu_pwv_radEXP xcor_pwv_radEXP TR_AAo_radEXP TR_AbdAo_radEXP ...
    ttp_pwv_radHR ttf_pwv_radHR ttu_pwv_radHR xcor_pwv_radHR TR_AAo_radHR TR_AbdAo_radHR};
basedir = 'I:\LIFE\patients\Visit1';
cd(basedir);

d=dir('*life*');
for i=9:length(d)
    disp(d(i).name)
    cd(d(i).name)
    life = findstr(d(i).name,'life');
    lifeid{i} = d(i).name(life(2):life(2)+8);
    cd('dicoms')
    dd=dir('*PWV_CartBH_AAo');
    if isempty(dd)
        scan_date{i}= 'UNK';
        cd(basedir)
    else
        cd(dd(1).name);
        ddd = dir('*dcm');
        if isempty(ddd)
            istgz = dir('*.tgz');
            if ~isempty(istgz)
                gunzip('*.tgz', 'Dicoms'); %unzip first
                cd('Dicoms') %move to unzipped folder
                dddd = dir(); %get the name of the only file in the new dir
                untar(dddd(3).name,'Dicoms'); %untar that file
                movefile('Dicoms/*','..'); %move unzipped files back up
                cd('..') %move up a directory
                rmdir('Dicoms','s') %get rid of created dummy unzipping folder
                ddd = dir('*dcm');
            end 
        end 
        info = dicominfo(ddd(1).name);
        dt = datetime([info.AcquisitionDate ' ' info.AcquisitionTime],'InputFormat','yyyyMMdd HHmmss');
        scan_date{i} = datestr(dt);
        mri_loc{i} = info.InstitutionName;
        coil{i} = info.ReceiveCoilName;

        lifeid = lifeid';
        scan_date = scan_date';
        mri_loc = mri_loc';
        coil = coil';
        cd(basedir)
    end 
        
%     a = dir('*radial_Aao');
%     if ~isempty(a)
%         cd(a.name);
%         aa = dir('10000p_LLR*');
%         if ~isempty(aa)
%             cd(aa.name);
%             cd('dat');
% 
%             fid = fopen('pcvipr_header.txt', 'r'); %open header
%             dataArray = textscan(fid,'%s%s%[^\n\r]','Delimiter',' ', ...
%                 'MultipleDelimsAsOne',true,'ReturnOnError',false); %parse header info
%             fclose(fid);
%             dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:),'UniformOutput',false);
%             pcviprHeader = cell2struct(dataArray{1,2}(:),dataArray{1,1}(:),1); %turn to structure
%             TR_AAo(i) = pcviprHeader.timeres;
%             cd(basedir);
%         end 
%     end 
%     cd(basedir)
%     
%     cd(d(i).name)
%     a = dir('*radial_AbdAo');
%     if ~isempty(a)
%         cd(a.name);
%         aa = dir('10000p_LLR*');
%         if ~isempty(aa)
%             cd(aa.name);
%             cd('dat');
% 
%             fid = fopen('pcvipr_header.txt', 'r'); %open header
%             dataArray = textscan(fid,'%s%s%[^\n\r]','Delimiter',' ', ...
%                 'MultipleDelimsAsOne',true,'ReturnOnError',false); %parse header info
%             fclose(fid);
%             dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:),'UniformOutput',false);
%             pcviprHeader = cell2struct(dataArray{1,2}(:),dataArray{1,1}(:),1); %turn to structure
%             TR_AbdAo(i) = pcviprHeader.timeres;
%             cd(basedir);
%         end   
%     end 
%     cd(basedir)
end 
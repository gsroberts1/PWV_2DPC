basedir = 'E:\LIFE\patients';
cd(basedir);

d=dir('*life*');
for i=1:length(d)
    disp(d(i).name)
    cd(d(i).name)
    life = findstr(d(i).name,'life');
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
        z1 = info.ImagePositionPatient(3);
    end 

    cd ..
    dd=dir('*PWV_CartBH_AbdAo');
    if ~isempty(dd)
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
        z2 = info.ImagePositionPatient(3);
    end 
    if (exist('z1','var') && exist('z2','var'))
        dz(i) = abs(z1-z2)';
        lifeid{i} = d(i).name(life(2):life(2)+8);
    end 
    clear z1 z2
    cd(basedir)
end 
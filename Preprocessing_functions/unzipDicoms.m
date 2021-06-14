% cd('F:\PWV');
% lifeDir = dir('*life*'); %get list of all life patient directories
% 
% for i=1:length(lifeDir) 
%     cd(lifeDir(i).name); %go into patient directory
%     cd('dicoms') %go into dicoms folder
%     dicomDir = dir(); %list all dicom folders
%     for j=3:length(dicomDir)
%         cd(dicomDir(j).name); %go into dicom folder
%         d = dir(); %list all files in dicom folder
%         if sum(contains({d.name},'dcm'))==0 %if we haven't unzipped...
%             gunzip('*.tgz', 'Dicoms'); %unzip first
%             cd('Dicoms') %move to unzipped folder
%             dd = dir(); %get the name of the only file in the new dir
%             untar(dd(3).name,'Dicoms'); %untar that file
%             movefile('Dicoms/*','..'); %move unzipped files back up
%             cd('..') %move up a directory
%             rmdir('Dicoms','s') %get rid of created dummy unzipping folder
%         end 
%         cd('..') %move back up to directory with all dicoms
%     end 
%     cd('../..') %move back to list of all patients
% end 


unzipDir = uigetdir();
cd(unzipDir)
dicomDir = dir(unzipDir); %get list of all life patient directories

    for j=3:length(dicomDir)
        cd(dicomDir(j).name); %go into dicom folder
        d = dir(); %list all files in dicom folder
        if sum(contains({d.name},'dcm'))==0 %if we haven't unzipped...
            gunzip('*.tgz', 'Dicoms'); %unzip first
            cd('Dicoms') %move to unzipped folder
            dd = dir(); %get the name of the only file in the new dir
            untar(dd(3).name,'Dicoms'); %untar that file
            movefile('Dicoms/*','..'); %move unzipped files back up
            cd('..') %move up a directory
            rmdir('Dicoms','s') %get rid of created dummy unzipping folder
        end 
        cd('..') %move back up to directory with all dicoms
    end 
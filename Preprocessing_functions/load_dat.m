function v = load_dat(name, res)
%% Load Dat
% Loads in dat files in current directory.
[fid,errmsg]= fopen(name,'r');
if fid < 0  % If name does not exist in directory
    disp(['Error Opening Data : ',errmsg]);
end

% Reads in as short, reshapes by image res.
temp = fread(fid,'short=>single');
v = reshape(temp,res(1),res(2));
fclose(fid);
end 
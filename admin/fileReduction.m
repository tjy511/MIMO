% cd to file
startup
fileLoc = strcat(rwd,'data/process/mimo/unattended/');
%cd(strcat(rwd,'data/process/mimo/unattended/'))
cd(fileLoc);

newLoc = strcat(fileLoc,'select/');

% List files to load
dirList = dir('*.mat');
for ii = 1:length(dirList)
    fileList{ii,1} = dirList(ii).name;
end

%% 
for ii = 1:length(fileList)
load(strcat(fileLoc,fileList{ii}))
save(strcat(newLoc,fileList{ii}),'xxPix','yyPix','imgPlane','vdat','Rs')
end
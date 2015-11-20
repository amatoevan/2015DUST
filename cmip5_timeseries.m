%% Calculate and save trends from projection of EOF #2 onto CMIP5 models
% Here this is for rcp8.5, can be for any run though

% Read in the CMIP5 data (from cmip5EOFproject.m)
system('ls rcp85/PC*rcp85*mat > rcp.list');
fileID=fopen('rcp.list');
files = textscan(fileID,'%s');
fclose(fileID);
files = files{1,1};

% load scale factors to conver the PC time series into units of DAOD
load('../DATA/scale_factors.mat')
nf = length(files);
y = NaN(1140/12,nf);
for i = 1:nf
    load(files{i})
    pc2 = mean(reshape(pc2,12,1140/12))';
    pc2 = (pc2-mean(pc2))/std(pc2);
    pc2 = pc2*sf + off;
    y(:,i) = pc2;
end

% save the time series
data.x = 2006+(1:1140/12)';
data.y = y;
names = files;
for i = 1:nf
    tmp = files{i};
    tmp = tmp(17:end-20);
    names{i} = tmp;
end
data.names = names;
save('out.rcp85.mat','-struct','data');

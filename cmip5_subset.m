% Read in and subset 10m winds from CMIP 5 simulations
% Only over N. Africa
% Make ensemble means if there is more than one simulation

% time of the rcp runs
time = 2006:1/12:(2100+11/12);

% list all model paths (for uas only)
system('ls uas/ > uas.list');
system('ls vas/ > vas.list');

% get file paths+names
fileID=fopen('uas.list');
m_pth=textscan(fileID,'%s');
fclose(fileID);
m_pth=m_pth{1,1};

% read in and make ensemble means
nmodel = length(m_pth);
for i = 1:nmodel
  pth_u = ['uas/' m_pth{i} '/'];
  system(['ls ' pth_u '*.nc > tmp.list']);
  fileID=fopen('tmp.list');
  fns_u=textscan(fileID,'%s');
  fclose(fileID);
  fns_u=fns_u{1,1};
  
  pth_v = ['vas/' m_pth{i} '/'];
  system(['ls ' pth_v '*.nc > tmp.list']);
  fileID=fopen('tmp.list');
  fns_v=textscan(fileID,'%s');
  fclose(fileID);
  fns_v=fns_v{1,1};
  
  nfls = length(fns_u);
  if nfls > 0

  for j = 1:nfls
    file_u = [fns_u{j}];
    file_v = [fns_v{j}];
    if j==1;
      y = ncread(file_u,'lat');
      x = ncread(file_u,'lon');
      u10 = NaN(length(x),length(y),length(time),nfls);
      v10 = NaN(length(x),length(y),length(time),nfls);
    end
    u10(:,:,:,j) = ncread(file_u,'uas');
    v10(:,:,:,j) = ncread(file_v,'vas');
  end
  u10 = nanmean(u10,4);
  v10 = nanmean(v10,4);

  % swap the longitudes
  ix1 = find(x<=180);
  ix2 = find(x>180);
  u10 = u10([ix2' ix1'],:,:);
  v10 = v10([ix2' ix1'],:,:);
  x = x([ix2' ix1']);
  x(x>180) = x(x>180)-360;

  % subset over N. Africa
  ix = x>=-20 & x<=40;
  iy = y>=0 & y<=40;
  u10 = u10(ix,iy,:);
  v10 = v10(ix,iy,:);
  x = x(ix);
  y = y(iy);

  % save the file (here for rcp 8.5)
  data.u10 = u10;
  data.v10 = v10;
  data.time = time;
  data.lon = x;
  data.lat = y;
  svnm = ['rcp85-' m_pth{i} '.nafrica.monthly.mat'];
  save(svnm,'-struct','data');

  end
end



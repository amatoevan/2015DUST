%% Project 10m wind EOF #2 from ERA-I data onto CMIP models

% Read in the subsetted cmip5 runs (from cmip5_subset.m)
system('ls rcp85*mat > rcp.list');
fileID=fopen('rcp.list');
files = textscan(fileID,'%s');
fclose(fileID);
files = files{1,1};
nf = length(files);

for ij = 1:nf
  % Load model output
  display(files{ij});
  load(files{ij});
  s = size(u10);
  u10 = double(u10(:,:,1:floor(s(3)/12)*12));
  v10 = double(v10(:,:,1:floor(s(3)/12)*12));
  time = (2006:1/12:(2100+11/12))';

  % spatially smooth the data (3x3)
  s = size(u10);
  u = u10;
  v = v10;
  for i = 1:s(3);
      u(:,:,i) = filter2(ones(3,3)/9,u(:,:,i));
      v(:,:,i) = filter2(ones(3,3)/9,v(:,:,i));
  end

  % remove the seasonal cycle
  u = reshape(u,s(1)*s(2),12,s(3)/12);
  v = reshape(v,s(1)*s(2),12,s(3)/12);
  for i = 1:s(1)*s(2);
   u(i,:,:)=squeeze(u(i,:,:))-repmat(nanmean(squeeze(u(i,:,:)),2),1,s(3)/12);
   v(i,:,:)=squeeze(v(i,:,:))-repmat(nanmean(squeeze(v(i,:,:)),2),1,s(3)/12);
  end
  u = reshape(u,s(1)*s(2),s(3));
  v = reshape(v,s(1)*s(2),s(3));

  % Change the resolution to match ERA-I
  % Load up the ERA-I EOF/PC results (identical to that calculated in
  % Figure1.m)
  EOF = load('../DATA/EOF.erai.mat');
  u = reshape(u,s(1),s(2),s(3));
  v = reshape(v,s(1),s(2),s(3));

  [Ex,Ey] = meshgrid(EOF.lon,EOF.lat);
  Ex = double(Ex)'; Ey = double(Ey)';

  [x,y] = meshgrid(lon,lat);
  x = double(x)'; y = double(y)';

  Uq = NaN(length(EOF.lon)*length(EOF.lat),s(3));
  Vq = NaN(length(EOF.lon)*length(EOF.lat),s(3));
  for i = 1:s(3);
      tmp = u(:,:,i);
      F = scatteredInterpolant(x(:),y(:),tmp(:));
      Uq(:,i) = F(Ex(:),Ey(:));
      tmp = v(:,:,i);
      F = scatteredInterpolant(x(:),y(:),tmp(:));
      Vq(:,i) = F(Ex(:),Ey(:));
  end

  s = size(Uq);
  Ex = Ex(:);
  Ey = Ey(:);

  % only N. Africa and over land
  Uq = Uq(EOF.land==1 & Ex<30 & Ey>5 & Ey<35,:);
  Vq = Vq(EOF.land==1 & Ex<30 & Ey>5 & Ey<35,:);
  uv = [Uq' Vq']';

  % project the cmip5 data onto the EOF pattern
  PCs = uv'*EOF.EOF;
  
  % save the output
  data.time = time;
  data.pc2 = PCs(:,2);
  save(['PC2.' files{ij}],'-struct','data')
  clearvars data

end


%% EOF analysis of 925hPa winds over N. Africa
% Should be:
% 1. land only: 5-35N & 20W-30E
% 2. Use monthly and de-seasonalized data
coast = load('coast.mat');
set(0,'defaulttextfontsize',6);


%% Remove seasonal cycle, get lat weights
load('~/Documents/Bodele/Modes/DATA/erai.nafrica.monthly.mat')
s = size(u);

% smooth the data
for i = 1:s(3);
    u(:,:,i) = filter2(ones(3,3)/9,u10(:,:,i));
    v(:,:,i) = filter2(ones(3,3)/9,v10(:,:,i));
    t(:,:,i) = filter2(ones(3,3)/9,t(:,:,i));
end

% remove the seasonal cycle
land = reshape(land,s(1)*s(2),1);
u = reshape(u,s(1)*s(2),12,s(3)/12);
v = reshape(v,s(1)*s(2),12,s(3)/12);
t = reshape(t,s(1)*s(2),12,s(3)/12);
for i = 1:s(1)*s(2);
 u(i,:,:)=squeeze(u(i,:,:))-repmat(nanmean(squeeze(u(i,:,:)),2),1,s(3)/12);
 v(i,:,:)=squeeze(v(i,:,:))-repmat(nanmean(squeeze(v(i,:,:)),2),1,s(3)/12);
 t(i,:,:)=squeeze(t(i,:,:))-repmat(nanmean(squeeze(t(i,:,:)),2),1,s(3)/12);
end
u = reshape(u,s(1)*s(2),s(3));
v = reshape(v,s(1)*s(2),s(3));
t = reshape(t,s(1)*s(2),s(3));

[x_,y_] = meshgrid(lon,lat);
x_ = x_'; y_ = y_';
x_ = reshape(x_,s(1)*s(2),1);
y_ = reshape(y_,s(1)*s(2),1);

cosW = repmat(cos(double(y_)),1,s(3));

u = u(land==1 & x_<30 & y_>5 & y_<35,:);
v = v(land==1 & x_<30 & y_>5 & y_<35,:);
t = t(land==1 & x_<30 & y_>5 & y_<35,:);
uv = [u' v']';
cosW = cosW(land==1 & x_<30 & y_>5 & y_<35,:);
len = length(land(land==1 & x_<30 & y_>5 & y_<35));

stduv = repmat(std(uv,0,2),1,s(3));
stdt = repmat(std(t,0,2),1,s(3));


%% DEM
lat75 = -flipud(lat75);

f1 = figure(1);
f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [0 0 8.9 8.9*.7];
f1.PaperPosition = [4.7 9.24 8.9 8.9*.7];
f1.PaperSize = [18.3 24.7];

cm = (colormap('bone'));
cm = cm(fliplr(1:64),:);
colormap(cm);

z = elv75;
[~,ch]=contourf(lon75,lat75,z',0:200:2000,'linewidth',.5);
set(gca,'ydir','normal','fontsize',6)
set(gca,'xlim',[-20 30],'ylim',[5 38])

cb = colorbar;
ylabel(cb,'Elevation (m)')
cb.FontSize = 6;
cb.Ticks = 0:250:2000;
for i = 1:length(cb.TickLabels);
    if length(cb.TickLabels{i}) > 3;
        cb.TickLabels{i} = [cb.TickLabels{i}(1) ',' cb.TickLabels{i}(2:4)];
    end
end

hold on
plot(coast.long,coast.lat,'-k','linewidth',1);
for i = 1:4;
    pp = plot([-20 -20]+10*i,[5 38],'color',[.5 .5 .5],'linewidth',.1);
    pp.Color(4)=.3;
end
for i = 1:6;
    pp = plot([-20 30],[5 5]+5*i,'color',[.5 .5 .5],'linewidth',.1);
    pp.Color(4)=.3; 
end

% Data from Evan et al. 2015 (Aeolian Research).
% NETCDF file available online at evan.ucsd.edu
f = '~/Documents/BODELE/DATA/seviri.emission.nc';
sev.lon = ncread(f,'lon');
sev.lat = ncread(f,'lat');
sev.emi = ncread(f,'emi');

s = size(sev.emi);
for i = 2:s(1)-1;
    for j = 2:s(2)-1;
        if sev.emi(i,j) >= .1;
            x = sev.lon(i)-.5;
            y = sev.lat(j)-.5;
            plot([x x+1],[y y+1],'b')
            if sev.emi(i,j+1) < .1 
                plot([x x+1],[y y],'b');
            end
            if sev.emi(i,j-1) < .1
                plot([x x+1],[y+1 y+1],'b');
            end
            if sev.emi(i-1,j) < .1
                plot([x x],[y y+1],'b');
            end
            if sev.emi(i+1,j) < .1
                plot([x+1 x+1],[y y+1],'b');
            end
        end
    end
end
xlabel('Longitude'); ylabel('Latitude');


% Elevation Transects
z = elv75;
z(z<0) = NaN;
[x,y] = meshgrid(lon75,lat75);x=x';y=y';

xT = NaN(10,4);
yT = NaN(10,4);
dT = NaN(10,4);
transT = NaN(10,4);

% Atlas Elevation transect
xo = -6.75;
yo = 31.5;
dx = 6.75+5.25;
dy = 23.25-31.5;
dr = .75;
theta = atan(abs(dx/dy));
dp = dr*sin(theta);
dq = dr*cos(theta);
hold on
fill([xo xo+dx xo+dx+dq xo+dq xo],[yo yo+dy yo+dy+dp yo+dp yo],'w')

% Air Elevation transect
xo = 5.25;
yo = 23.25;
dx = 17.25-5.25;
dy = 21-23.25;
theta = atan(abs(dx/dy));
dp = dr*sin(theta);
dq = dr*cos(theta);
fill([xo xo+dx xo+dx+dq xo+dq xo],[yo yo+dy yo+dy+dp yo+dp yo],'w')

% Hogar Air
xo = 5.25;
yo = 23.25;
dx = 9-5.25;
dy = 18.75-23.25;
theta = atan(abs(dx/dy));
dp = dr*sin(theta);
dq = dr*cos(theta);
fill([xo xo+dx xo+dx+dq xo+dq xo],[yo yo+dy yo+dy+dp yo+dp yo],'w')

% Bodele Elevation transect
xo = 17.25;
yo = 21;
dx = 24-17.25;
dy = 13.5-21;
theta = atan(abs(dx/dy));
dp = dr*sin(theta);
dq = dr*cos(theta);
fill([xo xo+dx xo+dx+dq xo+dq xo],[yo yo+dy yo+dy+dp yo+dp yo],'w')

F = scatteredInterpolant(x(:),y(:),z(:));
for i = 1:4;
    dT(:,i) = ...
        linspace(1,sqrt((max(xT(:,i))-min(xT(:,i))).^2 + ...
        (max(yT(:,i))-min(yT(:,i))).^2),10);
    transT(:,i) = F(xT(:,i),yT(:,i));
end


% Label the transects
labels = {'a','b','c','d','e'};
xT(1:5) = [-9.5 6 18.5 10.5 24];
yT(1:5) = [33 25.5 23 17.5 12];
for i = 1:5;
    t1 = text(xT(i),yT(i),labels(i));
    set(t1,'fontsize',7,'EdgeColor',[0 0 0],'BackgroundColor',[1 1 1]);
    set(t1,'VerticalAlignment','Middle');
end

leg = {'a) Atlas Mountains','b) Ahaggar Mountains', ...
       'c) Tibesti Mountains','d) Aïr Massif','e) Ennedi Plateau'};
tl = text(-18.3,9.5,leg);
set(tl,'fontsize',6,'EdgeColor',[0 0 0],'BackgroundColor',[1 1 1]);
set(tl,'VerticalAlignment','Middle');

hold off


%% print and write csv file
print(f1,'-dpdf','-cmyk','figure2.pdf')

[x,y] = meshgrid(lon75,lat75);
A = [reshape(x',80*54,1) reshape(y',80*54,1) reshape(z,80*54,1)];
csvwrite('figure2a1.csv',A)

[x,y] = meshgrid(sev.lon,sev.lat);
A = [reshape(x',65*35,1) reshape(y',65*35,1) reshape(sev.emi,65*35,1)];
csvwrite('figure2a2.csv',A)

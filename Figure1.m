%% Make Figure 1 of the Nature paper:
%  EOF analysis of 925hPa winds over N. Africa and plot along with the time
%  series of dust over Cape Verde.
%  Should be:
%    1. land only: 5-35N & 20W-30E
%    2. Use monthly and de-seasonalized data
%  Can be used for any reanalysis data set (modify line 12)
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


%% EOF: Do an EOF analysis and remove higher-order terms
C = (uv.*[cosW' cosW']') * (uv.*[cosW' cosW']')';
C = C/(s(3)-1);
[EOF,lam] = eigs(C,15); % saving 15 terms

sig = diag(lam);
B = uv'*EOF;

ruv = NaN(len*2,3);
for j = 1:3;
    x = B(:,j); x = x/std(x);
    for i = 1:len*2;
        b = regress(uv(i,:)',[ones(size(x)) x]);
        ruv(i,j) = b(2);
    end
end

% For the regression map
rt = NaN(len,1);
x = B(:,2); x = x/std(x);
for i = 1:len;
    b = regress(t(i,:)',[ones(size(x)) x]);
    rt(i,1) = b(2);
end

% split up the u & v components
ru = ruv(1:len,:);
rv = ruv((len+1):end,:);

% add back in the over-water NaN parts
tmp = NaN(s(1)*s(2),3);
tmp(land==1 & x_<30 & y_>5 & y_<35,:) = ru;
ru = reshape(tmp,s(1),s(2),3);

tmp = NaN(s(1)*s(2),3);
tmp(land==1 & x_<30 & y_>5 & y_<35,:) = rv;
rv = reshape(tmp,s(1),s(2),3);

tmp = NaN(s(1)*s(2),1);
tmp(land==1 & x_<30 & y_>5 & y_<35,1) = rt;
rt = reshape(tmp,s(1),s(2));


%% Plot SVD output
x_ = reshape(x_,s(1),s(2));
y_ = reshape(y_,s(1),s(2));

cm = colormap('parula');
cm(1,:) = ones(1,3);
colormap(cm);

f1 = figure(1);
f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [10 10 8.9 8.9*1.4];
f1.PaperPosition = [4.7 6.12 8.9 8.9*1.4];
f1.PaperSize = [18.3 24.7];

% Panel A
subplot(2,1,1)
z = sqrt(ru(:,:,2).^2+rv(:,:,2).^2);
z(isnan(z)) = 0;
imagesc(lon,lat,z')
cb = colorbar;
set(cb,'fontsize',6)
caxis([0 floor(max(z(:)))]);
grid on; box on
hold on
plot(coast.long,coast.lat,'k','linewidth',2);
caxis([0 floor(8*max(z(:)))/10 ])

% ru(1,:,2) = .5; rv(1,:,2) = 0;
q = quiver(x_(1:2:end,1:2:end),y_(1:2:end,1:2:end), ...
    ru(1:2:end,1:2:end,2),rv(1:2:end,1:2:end,2),2,'k');
hold off

% Text and position
set(gca,'xlim',[-20 30],'ylim',[5 35]);
set(gca,'ydir','normal','fontsize',6);
ylabel(cb,'Windspeed (ms^-^1 per std. dev.)')
title(['EOF/PC ' num2str(2) ': ' ...
    num2str(round(100*sig(2)/sum(sig))) '% Variance Exp.'], ...
    'FontWeight','Normal','fontsize',6);
colormap(cm);
xlabel('Longitude'); ylabel('Latitude');
text(-30.,35,'a','FontWeight','Bold','fontsize',8)
% set(gca,'position',[.08 .15 .3 .72])

% Legend for arrows
% axes('Position',[.34 0.03 .07 .09]);
axes('Position',[.7 .515 .2 .05]+[.07 0 -.105 0]);
q2 = quiver(1,1,1,0,1,'k');
set(gca,'xlim',[1.7 2.15],'ylim',[.93 1.07])
txt = text(.72,.95,'0.5 ms^{-1} per std. dev.');
axis off


% save wind speed & lat/lon values for .csv file
lon = x_; lat = y_;
uwind = ru(:,:,2); vwind = rv(:,:,2);


% Panel b
time = (1979:1/12:(2014+11/12))';
if s(3) < 200;
    time = 1979:1/4:(2014+3/4);
end
avhrr = load('~/Documents/Bodele/Modes/DATA/capeverde.avhrr.daod.mat');

subplot(2,1,2)
pc = smooth(B(:,2)./std(B(:,2)),13);
pc = (pc-mean(pc))/std(pc);
plot(time(7:end-6),pc(7:end-6)) % uv

hold on
nx = avhrr.time;
ny = smooth(avhrr.daod,13);
ny = (ny-mean(ny))/std(ny);
plot(nx(7:end-6),ny(7:end-6),'linewidth',1)

grid on
box on
set(gca,'xlim',[1979 2015],'ylim',[-2 3],'fontsize',6);
xlabel('Year');
ylabel('Standard deviation');
lh=legend('10m Winds','AVHRR dust');
set(lh,'location','northeast','fontsize',6)
plot(time,time*0,'k')
hold off

[r,p] = corr2(ny(7:end-6),pc(time>=nx(7)&time<=nx(end-6)));
ttl = title(['r(p)-value: ' num2str(round(r*100)/100) '' ...
    '(' num2str(round(p*100)/100) ')'],'FontWeight','Normal', ...
    'fontsize',6);
text(1974,3,'b','FontWeight','Bold','fontsize',8)
% set(gca,'position',[.63 .15 .32 .72])


%% Print & write
print(f1,'-dpdf','-cmyk','figure1.pdf')

A = [reshape(lon,61*41,1) reshape(lat,61*41,1) ...
    reshape(uwind,61*41,1) reshape(vwind,61*41,1)];
csvwrite('figure1a.csv',A)

A = [avhrr.time avhrr.daod];
csvwrite('figure1b1.csv',A)

A = [time B(:,2)./std(B(:,2))];
csvwrite('figure1b2.csv',A)

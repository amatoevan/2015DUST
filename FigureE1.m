%% Make Figure E1 of the Nature paper:
%  EOF analysis of 10m winds over N. Africa
%  Should be:
%    1. land only: 5-35N & 20W-30E
%    2. Use monthly and de-seasonalized data
%  Can be used for any reanalysis data set (modify line 12)
coast = load('coast.mat');
set(0,'defaulttextfontsize',6);


% Remove seasonal cycle, get lat weights
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


% EOF: Do an EOF analysis and remove higher-order terms
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
f1.Position = [0 0 18.3 10.];
f1.PaperPosition = [0 7.35 18.3 10.];
f1.PaperSize = [18.3 24.7];

ru_ = ru;
rv_ = rv;
% EOF Map
ind = {'a','b','c'};
for i = 1:3
    subplot(2,3,i)
    z = sqrt(ru(:,:,i).^2+rv(:,:,i).^2);
    z(isnan(z)) = 0;
    imagesc(lon,lat,z')
    cb = colorbar;
    set(cb,'fontsize',6)
    caxis([0 floor(max(z(:)))]);
    grid on; box on
    hold on
    plot(coast.long,coast.lat,'k','linewidth',2);
    
    di = 3;
    ru(:,:,i) = ru(:,:,i)./z;
    rv(:,:,i) = rv(:,:,i)./z;
    q = quiver(x_(1:di:end,1:di:end),y_(1:di:end,1:di:end), ...
        ru(1:di:end,1:di:end,i),rv(1:di:end,1:di:end,i),.8,'k');
    hold off
    
    % Text and position
    set(gca,'xlim',[-20 30],'ylim',[5 35]);
    set(gca,'ydir','normal','fontsize',6);
    ylabel(cb,'Windspeed (ms^-^1 per std. dev.)')
    ttl = title(['EOF/PC ' num2str(i) ': ' ...
        num2str(round(100*sig(i)/sum(sig))) '% Variance Exp.'] ,...
        'FontWeight','Normal','fontsize',6);
    ttl.Position = ttl.Position + [0 1 0];
    colormap(cm);
    xlabel('Longitude'); ylabel('Latitude');
    text(-37,35,ind(i),'FontWeight','Bold','fontsize',8)
    
    caxis([0 floor(8*max(z(:)))/10 ])
    
end


% Panel b
time = (1979:1/12:(2014+11/12))';
if s(3) < 200;
    time = 1979:1/4:(2014+3/4);
end

ind = {'d','e','f'};
for i = 1:3;
    subplot(2,3,i+3)
    pc = smooth(B(:,i)./std(B(:,i)),13);
    pc = (pc-mean(pc))/std(pc);
    plot(time(7:end-6),pc(7:end-6)) % uv
    
    hold on
    grid on
    box on
    set(gca,'xlim',[1979 2015],'ylim',[-2 3],'YTick',-2:3,'fontsize',6);
    xlabel('Year');
    ylabel('Standard deviation');
    plot(time,time*0,'k')
    hold off
    
    title(['PC ' num2str(i) ],'FontWeight','Normal','fontsize',6);
    text(1971,3,ind(i),'FontWeight','Bold','fontsize',8)
    set(gca,'position',get(gca,'position')+[0 0 -.05 0])
end


%% Print & write
print(f1,'-dpdf','-cmyk','figureE1.pdf')

[x,y] = meshgrid(lon,lat);
lon = x'; lat = y';
A = [reshape(lon,61*41,1) reshape(lat,61*41,1) ...
    reshape(ru_(:,:,1),61*41,1) reshape(rv_(:,:,1),61*41,1) ...
    reshape(ru_(:,:,2),61*41,1) reshape(rv_(:,:,2),61*41,1) ...
    reshape(ru_(:,:,3),61*41,1) reshape(rv_(:,:,3),61*41,1) ];
csvwrite('figureE1abc.csv',A)

A = [time B(:,1:3)];
csvwrite('figureE1def.csv',A)


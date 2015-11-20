%% Project EOF #2 onto Reanalysis and plot, along with other time series
% Then plot time series of CMIP5 EOF2 projections from RCP 8.5
% Figure 3 from Evan et al. 2015 Nature
set(0,'defaulttextfontsize',6);
coast = load('coast.mat');

% time span for standardization
rng = [1982 2008];

% AVHRR daod data (available at evan.ucsd.edu)
avhrr = load('~/Documents/Bodele/Modes/DATA/capeverde.avhrr.daod.mat');
ax = (1982:2009)';
ay = nanmean(reshape(avhrr.daod,12,28))';
ind = ax >= rng(1) & ax <= rng(2);
ai = ind;

% CIRES 20CR
load('~/Documents/Bodele/Modes/DATA/cires.pc2.mat');
s = size(pc2);
pc2 = reshape(pc2,12,s(1)/12);
pc2 = nanmean(pc2)';
cires.time = floor(min(time)):floor(max(time));
ref = cires.time>=rng(1)&cires.time<=rng(2);
cires.pc2 = (pc2-mean(pc2(ref)))/std(pc2(ref));
cires.pc2x = cires.pc2;
cires.pc2x(cires.time<1950) = cires.pc2x(cires.time<1950)-1.6392;
cires.pc2 = cires.pc2*std(ay(ai))+mean(ay(ai));
cires.sf = std(cires.pc2); % std dev for scaling CMIP5 models
cires.off = mean(cires.pc2(end-14:end)); % std dev for scaling CMIP5 models

% EOF ERAI
load('~/Documents/Bodele/Modes/DATA/erai.pc2.mat');
s = size(pc2);
pc2 = reshape(pc2,12,s(1)/12);
pc2 = nanmean(pc2)';
erai.time = floor(min(time)):floor(max(time));
ref = erai.time>=rng(1)&erai.time<=rng(2);
erai.pc2 = (pc2-mean(pc2(ref)))/std(pc2(ref));
erai.pc2 = erai.pc2*std(ay(ai))+mean(ay(ai));
erai.off = mean(erai.pc2(end-15:end));

% Read in the CV data; 15N & 23W (available at evan.ucsd.edu)
f = '~/Documents/COAUTHOR/2015_AB_etal/DATA/F1/tnatl.historical.dust.nc';
cy = ncread(f,'dust_optical_depth');
cy = squeeze(cy(15,42,:));
cy = mean(reshape(cy,12,54))';
cx = (1955:2008)';
ind = cx >= rng(1) & cx <= rng(2);

% Read in the Barbados data (available from Joe Prospero)
raw = dlmread('~/Documents/COAUTHOR/2015_AB_etal/DATA/S1/Barbados.txt');
bx  = raw(:,1); by = raw(:,2);
by(6) = (by(5)+by(7))/2;
by = (by-mean(by))/std(by);
by = by*std(cy)+mean(cy);


%% Annual time series with a 1-3-7-3-1 smoothing
f1 = figure(1);
% f1.Units = 'centimeters';
% f1.Position = [10 10 8.9 12];
% f1.PaperUnits = 'centimeters';
% f1.PaperPosition = [0 0 8.9 12];
% f1.PaperSize = [8.9 12];

f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [10 10 8.9 8.9*1.4];
f1.PaperPosition = [4.7 6.12 8.9 8.9*1.4];
f1.PaperSize = [18.3 24.7];

subplot(2,1,1)

% smooth
cires.pc2 = smooth(smooth(cires.pc2,3),3);
cires.pc2x = smooth(smooth(cires.pc2x,3),3);
cires.time = cires.time(3:end-2); cires.pc2 = cires.pc2(3:end-2);
cires.pc2x = cires.pc2x(3:end-2);

erai.pc2 = smooth(smooth(erai.pc2,3),3);
erai.time = erai.time(3:end-2); erai.pc2 = erai.pc2(3:end-2);

ay_ = ay; % unsmoothed, for CMIP5 stuff
ax_ = ax; % unsmoothed
ay = smooth(smooth(ay,3),3);
ax = ax(3:end-2); ay = ay(3:end-2);

cy = smooth(smooth(cy,3),3);
cx = cx(3:end-2); cy = cy(3:end-2);

by = smooth(smooth(by,3),3);
bx = bx(3:end-2); by = by(3:end-2);

pa = plot(cires.time,cires.pc2,erai.time,erai.pc2,ax,ay,cx,cy);
hold on
plot([0 1],[0 0],'Color',[.7 .7 .7]);
plot([0 1],[0 0],'k','linewidth',1);
lh = legend('PC2 20CR (this study)', ...
    'PC2 ERA-I (this study)',...
    'AVHRR (Ref. 2)',...
    'Cape Verde (Ref. 2)');
set(lh,'location','northwest','fontsize',6);
lh.Position = lh.Position + [-0.03 .02 0 0];

grid on
set(gca,'fontsize',6,'xlim',[min(cires.time) max(cires.time)+11/12])
xlabel('Years')
ylabel('Dust optical depth (unitless)')

set(gca,'xlim',[1850 2015],'xtick',1850:25:2015)
set(gca,'ylim',[0.19 0.65])

pn = plot(cires.time,cires.pc2,erai.time,erai.pc2,ax,ay,cx,cy);
for i = 1:4
    pn(i).Color = pa(i).Color;
end

plot([1850 2015],[1 1]*mean(cires.pc2),'--','color',get(pa(1),'Color'))
text(1830,0.65,'a','fontsize',8,'FontWeight','Bold')

% for the csv file
a.cires_x = cires.time';
a.cires_y = cires.pc2;
a.erai_x = erai.time';
a.erai_y = erai.pc2;
a.avhrr_x = ax;
a.avhrr_y = ay;
a.capeverde_x = cx;
a.capeverde_y = cy;


%% Process the CMIP5 data
cmip = load('out.rcp85.mat');

% First scale the ensemble members
s = size(cmip.y);
cmipy.ys = cmip.y;
for i = 1:s(2);
    cmip.y(:,i) = (cmip.y(:,i)-mean(cmip.y(:,i)))/std(cmip.y(:,i));
    cmip.y(:,i) = cmip.y(:,i)*cires.sf;
    cmip.y(:,i) = cmip.y(:,i) - mean(cmip.y(1:5,i)) + mean(cires.pc2);
    cmip.ys(:,i) = smooth(smooth(cmip.y(:,i),3),3);
end

% % Do ensemble averages
% [unames,~]=unique(cmip.names);
% y_ = NaN(95,length(unames));
% ys_ = NaN(95,length(unames));
% for i = 1:length(unames);
%     index = strncmp(unames(i),cmip.names,length(unames{i}));
%     y_(:,i) = mean(cmip.y(:,index),2);
%     ys_(:,i) = mean(cmip.ys(:,index),2);
% end
% cmip.y = y_;
% cmip.ys = ys_;

cmip.ys = cmip.ys(3:end-2,:);
cmip.xs = cmip.x(3:end-2);
% cmip.names = names;

subplot(2,1,2)

plot(cmip.xs,cmip.ys(:,1),'Color',[.7 .7 .7])
hold on
plot(cmip.xs,mean(cmip.ys,2),'k')
x = cmip.x; y = mean(cmip.y,2);
[b,bint] = regress(y,[ones(size(x)) x]);
plot(x,x*b(2)+b(1),'r')
lh = legend('CMIP5 Models (RCP 8.5)','CMIP5 Multimodel Mean', ...
    'CMIP5 trend');
set(lh,'location','northwest','fontsize',6);
lh.Position = lh.Position + [-0.03 .02 0 0];

plot(cmip.xs,cmip.ys,'Color',[.7 .7 .7])
plot(cmip.xs,mean(cmip.ys,2),'k','linewidth',1)
plot(x,x*b(2)+b(1),'r')

%display([num2str(b(2)*100) '+/-' num2str(mean(abs(bint(2,:)-b(2)))*100)]);

grid on
xlabel('Years')
ylabel('Dust optical depth (unitless)')

% The rest of the plot
set(gca,'xlim',[2005 2100],'xtick',2000:20:2100,'fontsize',6, ...
    'ylim',[0.19 0.65])
plot([2005 2100],[1 1]*mean(cires.pc2),'--','color',get(pa(1),'Color'))
text(1993.2,0.65,'b','fontsize',8,'FontWeight','Bold')

hold off
 

%% print and write csv file
print(f1,'-dpdf','-cmyk','figure3.pdf')

x = (min(a.cires_x):max(a.erai_x))';
y = NaN(length(x),4);
y(x<=max(a.cires_x),1) = a.cires_y;
y(x>=min(a.erai_x),2) = a.erai_y;
y(x>=min(a.avhrr_x)&x<=max(a.avhrr_x),3) = a.avhrr_y;
y(x>=min(a.capeverde_x)&x<=max(a.capeverde_x),3) = a.capeverde_y;
csvwrite('figure3a.csv',[x y])

csvwrite('figure3b.csv',[cmip.xs cmip.ys ])

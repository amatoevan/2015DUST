% Make a plot of the correlations of PC 2 and 
%  1. Sahel rainfall in the summer
%  2. Jones NAO in the winter
%  3. Owen's thing in the summer (and winter?)
%  4. Natalie's PDSI in the summer
%  5. Nino 3.4 in the winter

%% Read in the data
nyr = 30; % # of year - 1 for the correlation

% NAO - statistically negatively correlated
load('~/Documents/Bodele/Modes/DATA/cires.pc2.mat')
raw = reshape(pc2,12,1932/12);
pc2 = mean(reshape(pc2,12,1932/12))';
time = (min(time):max(floor(time)))';
nao = load('~/Documents/Bodele/Modes/DATA/jones.nao.mat');
nao.nao = mean(reshape(nao.nao,12,2304/12))';
nao.time = (min(nao.time):max(floor(nao.time)))';
g = nao.time>=min(time)&nao.time<=max(time);
% [r,p]=corr(pc2,nao.nao(g));


% Strongest pre-1930 and 1994-96
x = pc2; y = nao.nao(g); t = time;
r = NaN(length(x),1); p = r;
for i = (nyr/2+1):(length(x))-nyr/2;
    [r_,p_]=corr( x(i-nyr/2:i+nyr/2),y(i-nyr/2:i+nyr/2) );
    r(i) = r_; p(i) = p_;
end
rNAO = r; pNAO = p; tNAO = t+15;


% ENSO - not stationary, correlated with Nina in 60s/70s 
tmp = dlmread('~/Documents/Bodele/Modes/DATA/nino34/nino34.txt');
nino.nino = mean(tmp(:,2:end),2);
nino.time = tmp(:,1);
gP = time >= 1870;
gN = nino.time <= 2011;
% [r,p] = corr(pc2(gP),nino.nino(gN));

% correlated with Nina in 1959?1971
x = pc2(gP); y = nino.nino(gN); t = time(gP);
r = NaN(length(x),1); p = r;
for i = (nyr/2+1):(length(x))-nyr/2;
    [r_,p_]=corr( x(i-nyr/2:i+nyr/2),y(i-nyr/2:i+nyr/2) );
    r(i) = r_; p(i) = p_;
end
rENSO = r; pENSO = p; tENSO = t+15;


% PDSI
nat = load('~/Documents/Bodele/Modes/DATA/natalie/natalie.mat');
nat = load('~/Documents/Bodele/Modes/DATA/pdsi.mat');
g = time>=min(nat.time)&time<=max(nat.time);

x = pc2(g); y = nat.PDSI; t = time(g);
r = NaN(length(x),1); p = r;
for i = (nyr/2+1):(length(x))-nyr/2;
    [r_,p_]=corr( x(i-nyr/2:i+nyr/2),y(i-nyr/2:i+nyr/2) );
    r(i) = r_; p(i) = p_;
end
rPDSI = r; pPDSI = p; tPDSI = t+15;


% Owen
owen = load('~/Documents/Bodele/Modes/DATA/Owen/owen.lat.mat');
g = time>=min(owen.time)&time<=max(owen.time);
x = pc2(g); y = mean(owen.NCEP)'; t = time(g);
r = NaN(length(x),1); p = r;
for i = (nyr/2+1):(length(x))-nyr/2;
    [r_,p_]=corr( x(i-nyr/2:i+nyr/2),y(i-nyr/2:i+nyr/2) );
    r(i) = r_; p(i) = p_;
end
rOWEN = r; pOWEN = p; tOWEN = t+15;


% SHL
hl = ncread('~/Documents/Bodele/DATA/hl_day_lag.nc','val');
ltm = mean(mean(hl,2));
hl_s = mean(hl); % seasonal cycle
hl_s = smooth([hl_s hl_s hl_s],21,'moving');
hl_s = (hl_s(366:365+365))';
hl_a = hl - repmat(hl_s,34,1) + ltm;
hl = nanmean(hl_a(1:end-1,121:243),2);
% hlt = (1979:2011)';

g = time>=1979;
x = mean(raw(5:9,g))'; y = hl; t = time(g);
r = NaN(length(x),1); p = r;
for i = (nyr/2+1):(length(x))-nyr/2;
    [r_,p_]=corr( x(i-nyr/2:i+nyr/2),y(i-nyr/2:i+nyr/2) );
    r(i) = r_; p(i) = p_;
end
rSHL = r; pSHL = p; tSHL = t+15;


%% Plot
f1 = figure(1);
f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [0 0 12 7.5];
f1.PaperPosition = [3.15 8.6 12 7.5];
f1.PaperSize = [18.3 24.7];


pm = plot(tNAO,rNAO,tENSO,rENSO,tPDSI,rPDSI,tOWEN,rOWEN,tSHL,rSHL);
pm(5).Color = [0 0 0];
colors = NaN(5,3);
for i = 1:5
    colors(i,:) = pm(i).Color;
end


lh = legend('Jones NAO','Niño 3.4','Sahel PDSI', ...
    'COA Latitude','Z_{SHL}');
set(lh,'fontsize',6)

pv = .1;
hold on
scatter(tNAO(pNAO<pv),rNAO(pNAO<pv),25,colors(1,:),'fill')
scatter(tENSO(pENSO<pv),rENSO(pENSO<pv),25,colors(2,:),'fill')
scatter(tPDSI(pPDSI<pv),rPDSI(pPDSI<pv),25,colors(3,:),'fill')
scatter(tOWEN(pOWEN<pv),rOWEN(pOWEN<pv),25,colors(4,:),'fill')
scatter(tSHL(pSHL<pv),rSHL(pSHL<pv),25,colors(5,:),'fill')

hold off
grid
set(gca,'xtick',1860:10:2010,'fontsize',6,'xlim',[1880 2011])
ylabel(['r-value (prior ' num2str(nyr+1) '-years)'])

%% Print figure and data
print(f1,'-dpdf','-cmyk','figureE4.pdf')

X = (min(tNAO):max(tSHL))';
Y = NaN(length(X),5);
y(X>=min(tNAO)&X<=max(tNAO),2) = rNAO;
y(X>=min(tENSO)&X<=max(tENSO),2) = rENSO;
y(X>=min(tPDSI)&X<=max(tPDSI),2) = rPDSI;
y(X>=min(tOWEN)&X<=max(tOWEN),2) = rOWEN;
y(X>=min(tSHL)&X<=max(tSHL),2) = rSHL;


A = [X Y];
csvwrite('figureE4.csv',A)


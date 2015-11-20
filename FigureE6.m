% Compare trends in DMP and PC2 for RCP8.5 
dmp = load('~/Documents/BODELE/CMIP5/rcp85/dust.rcp85.atl.mat');
pc2 = load('out.rcp85.mat');

% get the trends
dmp.b = NaN(length(dmp.names),1);
dmp.err = NaN(length(dmp.names),1);
for i = 1:length(dmp.names) 
    y = mean(reshape(dmp.dmp(i,:),12,94))';
    y = (y-mean(y))/std(y);
    x = (1:94)';
    [b,bint]=regress(y,[ones(size(x)) x]);
    dmp.b(i) = b(2);
    dmp.err(i) = mean(abs(b(2)-bint(2,:)));
end

pc2.b = NaN(length(pc2.names),1);
pc2.err = NaN(length(pc2.names),1);
for i = 1:length(pc2.names) 
    y = pc2.y(:,i);
    y = (y-mean(y))/std(y);
    x = (1:95)';
    [b,bint]=regress(y,[ones(size(x)) x]);
    pc2.b(i) = b(2);
    pc2.err(i) = mean(abs(b(2)-bint(2,:)));
end

% link the models
x = dmp.b;
xe = dmp.err;
y = NaN(size(x));
ye = NaN(size(x));
for i = 1:length(x)
    ind = strcmp(dmp.names{i},pc2.names);
    if find(ind)
        y(i) = pc2.b(ind);
        ye(i) = pc2.err(ind);
    end
end
g = ~isnan(y);
x = x(g)*100; y = y(g)*100;

%% Make the figure
f1 = figure(1);
f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [0 0 12 12*.6];
f1.PaperPosition = [3.15 8.75 12 12*.6];
f1.PaperSize = [18.3 24.7];


cm = colormap('parula');
col = round(linspace(1,64,length(x)));
sym = ('o+*xsd^<>pho+*xsd');

sc = gscatter(x,y,1:length(x));
for i = 1:length(x);
    sc(i).Color = cm(col(i),:);
    sc(i).MarkerFaceColor = cm(col(i),:);
    sc(i).Marker = sym(i);
    sc(i).MarkerSize = 5;
end

lh = legend(dmp.names(g));
set(lh,'location','eastoutside','fontsize',6)

set(gca,'fontsize',6)
grid on; box on
hold on
b = regress(y,[ones(size(x)) x]);
x_ = [-3 4];
plot(x_,x_,'--k')
plot(x_,x_*b(2),'r')
hold off
xlabel('Dust mass path trend (std. dev. per 100 yrs)')
ylabel('PC2 (std. dev. per 100)')
title('RCP 8.5')


%% Print figure and data
print(f1,'-dpdf','-cmyk','figureE6.pdf')


A = [x y];
csvwrite('figureE6.csv',A)

names = dmp.names(g);
for i = 1:length(names);
    display(names{i})
end



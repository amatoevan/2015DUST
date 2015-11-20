% rMake bar charts for Figure E8 from Evan et al. 2015 Nature paper

cd '~/Documents/BODELE/Modes/Nature_r1/DATA/';
[r2,reanalysis,station] = importfile( ...
    'r2_ws10m_synop_reanalyses_2000-2013.txt');
cd '~/Documents/BODELE/Modes/Nature_r2/FIGS/'

colormap('gray')

f1 = figure(1);
f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [0 0 12 12*.6];
f1.PaperPosition = [3.15 8.75 12 12*.6];
f1.PaperSize = [18.3 24.7];

r2 = r2(1:end-4);
s = size(r2);
y = reshape(r2,4,s(1)/4)';
b = bar(y);

set(gca,'xticklabel',strtrim(station(1:4:s(1))));
set(gca,'fontsize',6,'xlim',[.5 s(1)/4+.5])
grid on
ylabel('r^2')

reanalysis(2) = {'NNRP'};
reanalysis(3) = {'NCEP-DOE'};
reanalysis(4) = {'ERAI'};
lh = legend(strtrim(reanalysis(1:4)));
set(lh,'fontsize',6,'location','northwest')

%% Print & write
print(f1,'-dpdf','-cmyk','figureE8.pdf')

csvwrite('figureE8.csv',y)
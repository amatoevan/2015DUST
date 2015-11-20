% Elevation Transects: Figure E3 of Nature paper
load('~/Documents/Bodele/Modes/DATA/erai.nafrica.monthly.mat')

f1 = figure(1);
f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [0 0 12 12];
f1.PaperPosition = [3.15 6.35 12 12];
f1.PaperSize = [18.3 24.7];

lat75 = -flipud(lat75);
z = elv75;
z(z<0) = NaN;
[x,y] = meshgrid(lon75,lat75);x=x';y=y';

xT = NaN(10,4);
yT = NaN(10,4);
dT = NaN(10,4);
transT = NaN(10,4);

% Atlas Elevation transect
xT(:,1) = linspace(-6.75,5.25,10);
yT(:,1) = linspace(31.5,23.25,10);

% Air Elevation transect
xT(:,2) = linspace(5.25,17.25,10);
yT(:,2) = linspace(23.25,21,10);

% Hogar Air
xT(:,3) = linspace(5.25,9,10);
yT(:,3) = linspace(23.25,18.75,10);

% Bodele Elevation transect
xT(:,4) = linspace(17.25,24,10);
yT(:,4) = linspace(21,13.5,10);

F = scatteredInterpolant(x(:),y(:),z(:));
for i = 1:4;
    dT(:,i) = ...
        linspace(1,sqrt((max(xT(:,i))-min(xT(:,i))).^2 + ...
        (max(yT(:,i))-min(yT(:,i))).^2),10);
    transT(:,i) = F(xT(:,i),yT(:,i));
end

% Get matlab colors
% pp = plot(transT,'-','linewidth',2);
% col = [get(pp(1),'Color')' get(pp(2),'Color')' ...
%     get(pp(3),'Color')' get(pp(4),'Color')'];

% Bar plot the transects
ttls = {'Atlas to Ahaggar','Ahaggar to Tibesti', ...
    'Ahaggar to Aïr','Tibesti to Ennedi'};
ind = {'a','b','c','d'};
for i = 1:4
    subplot(2,2,i)
    bar(111.12*dT(:,i),transT(:,i),'FaceColor',[.7 .7 .7])
    grid on
    set(gca,'fontsize',6)
    set(gca,'xlim',111.12*[min(dT(:,i))-0 max(dT(:,i))+0]);
    xlabel('Distance (km)')
    ylabel('Elevation (m)')
    title(ttls(i),'fontsize',6,'fontweight','normal')
    txt = text(-.26,1,ind(i),'fontsize',8,'fontweight','bold', ...
        'units','Normalized');
    
    ytl = get(gca,'YTickLabel');
    for j = 1:length(ytl);
        if length(ytl{j}) > 3;
            ytl{j} = [ytl{j}(1) ',' ytl{j}(2:4)];
        end
    end
    set(gca,'YTickLabel',ytl)

    xtl = get(gca,'XTickLabel');
    for j = 1:length(xtl);
        if length(xtl{j}) > 3;
            xtl{j} = [xtl{j}(1) ',' xtl{j}(2:4)];
        end
    end
    set(gca,'XTickLabel',xtl)

end

% print and write csv file
print(f1,'-dpdf','-cmyk','figureE3.pdf')

A = [dT(:,1) transT(:,1) dT(:,2) transT(:,2) ...
    dT(:,3) transT(:,3) dT(:,4) transT(:,4)];
csvwrite('figureE3abcd.csv',A)

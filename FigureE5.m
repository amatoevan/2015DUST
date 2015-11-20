%% Plot trends from projection of EOF #2 onto CMIP5 models (RCP 8.5 & 4.5)
% 1. land only: 5-35N & 20W-30E
% 2. Use monthly and de-seasonalized data

% Make the Figure
f1 = figure(1);
f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [10 10 18.3 18.3];
f1.PaperPosition = [0 3.2 18.3 18.3];
f1.PaperSize = [18.3 24.7];

ind = {'a','b'};
for i = 1:2
    load(['out.rcp' num2str(8-4*(i-1)) '5.mat'])
    
    nf = length(names);
    b = NaN(nf,1);
    bint = NaN(2,1);
    err = NaN(2,1);
    for j = 1:nf
        [c,cint] = regress(y(:,j),[ones(size(x)) x/100]);
        b(j) = c(2);
        bint(:,j) = cint(2,:);
        err(j) = mean(abs(c(2)-cint(2,:)));
    end
    [c,cint] = regress(mean(y,2),[ones(size(x)) x/100]);
    c = c(2);
    cint = cint(2,:);
    cerr = mean(abs(c-cint));
    
    subplot(2,1,i)
    errorbar([c b']',[cerr err']','o')
    hold on
    plot([0 nf+2],[0 0],'-r')
    plot([0 nf+2],[0 0]+c,'--b')
    hold off
    grid on
    
    set(gca,'xlim',[0.5 nf+1.5],'xtick',1:(nf+1),'fontsize',6)
    title(['RCP ' num2str(8-4*(i-1)) '.5']);
    ylabel('\DeltaDust Optical Depth (100 yrs^-^1)');
    
    xlab{1} = 'Multimodel Mean';
    for j = 1:nf;
        xlab{j+1} = names(j);
    end
    hText = xticklabel_rotate(1:(nf+1),90,xlab);
    for j = 1:nf+1;
        hText(j).FontSize = 6;
        hText(j).Position = hText(j).Position+[-.005 -.03 0];
    end
    
    yl = get(gca,'YLim');
    text(-1.75,yl(2),ind(i),'FontWeight','Bold','fontsize',8)
    clearvars xlab
    
end

%% Print images and data
print(f1,'-dpdf','-cmyk','figureE5.pdf')

load('out.rcp45.mat')
A = [x y];
csvwrite('figureE5a.csv',A)
for i = 1:length(names);
    display(names{i})
end

display('----')

load('out.rcp85.mat')
A = [x y];
csvwrite('figureE5b.csv',A)
for i = 1:length(names);
    display(names{i})
end


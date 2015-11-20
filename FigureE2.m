%% Compare monthly mean wspd to DUP
%   P         4-D             2571264  double              
%   T         4-D             2571264  double              
%   dir       4-D             2571264  double              
%   mr        4-D             2571264  double              
%   rh        4-D             2571264  double              
%   wsp       4-D             2571264  double              
location = {'DAAJ','DAAT','DRZA'};
ttl = {'Djanet, Algeria','Tamanrasset, Algeria','Agadez, Niger'};


%% Scatter plots
f1 = figure(1);
f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [0 0 18.3 4.27];
f1.PaperPosition = [0 10.2 18.3 4.27];
f1.PaperSize = [18.3 24.7];

ind = {'a','b','c'};
for k = 1:3;
    subplot(1,3,k)
    
    load(['~/Documents/BODELE/Modes/Nature_r1/DATA/' location{k} '.mat'])
    U = wsp;
    U(U==0) = NaN;
    U(U>20) = NaN;
    E = U;
    
    %E = E.^3;
    E = (E.^3).*(1+6./E).*(1-36./(E.^2));
    E(U<6) = 0;
    
    U = nanmean(U,4);
    E = nanmean(E,4);
    
    s = size(U);
    for i = 1:s(1);
        for j = 1:s(2)     
            if length(find(~isnan(U(i,j,:)))) < 15;
                U(i,j,:) = NaN;
            end
        end
    end
    E(isnan(U)) = NaN;
    
    U = reshape(nanmean(U,3)',s(1)*s(2),1);
    E = reshape(nanmean(E,3)',s(1)*s(2),1);
    
    scatter(U,E,'fill')
    [b,bint]=regress(E,[ones(size(U)) U]);
    xl = get(gca,'xlim');
    hold on
    x = xl(1):.1:xl(2);
    plot(x,x*b(2)+b(1))
    hold off
    grid
    box on
    set(gca,'fontsize',6)
    xlabel('Wind speed (m s^-^1)')
    ylabel('DUP (unitless)')
    title(ttl{k},'fontsize',6,'fontweight','normal')

    set(gca,'ylim',[0 max(E)]);
    set(gca,'xlim',[min(U) max(U)]);
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    
    if k==3;
        set(gca,'YTick',0:200:1000);
        ytl = get(gca,'YTickLabel');
        ytl{end} = '1,000';
        set(gca,'YTickLabel',ytl)
    end

    [r,p] = corr(U(~isnan(U)),E(~isnan(U)));
    t1 = text( (xl(2)-xl(1))*.1+xl(1), (yl(2)-yl(1))*.9+yl(1), ...
        ['r = ' num2str( round(r*100)/100 ) ', p < 0.01']);
    t1.FontSize = 6;
    t1.BackgroundColor = [1 1 1];
    t1.EdgeColor = [0 0 0];
    
    text(xl(1)-(xl(2)-xl(1))*.22,yl(2),ind(k), ...
        'FontWeight','Bold','fontsize',8)
    
    if k==1;
        A.U = NaN(432,3);
        A.E = NaN(432,3);
    end
    A.U(:,k) = U;
    A.E(:,k) = E;
end


%% Print & write
print(f1,'-dpdf','-cmyk','figureE2.pdf')

A = [A.U(:,1) A.E(:,1) A.U(:,2) A.E(:,2) A.U(:,3) A.E(:,3)];
csvwrite('figureE2abc.csv',A)


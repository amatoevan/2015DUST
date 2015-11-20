%% Make Figure E7 of the Nature paper:
%  EOF analysis of 10m winds over N. Africa
%  Should be:
%    1. land only: 5-35N & 20W-30E
%    2. Use monthly and de-seasonalized data
%  Can be used for any reanalysis data set (modify line 12)
coast = load('coast.mat');
set(0,'defaulttextfontsize',6);

files = {'~/Documents/BODELE/Modes/DATA/EOF_CIRES.mat', ...
    '~/Documents/BODELE/Modes/DATA/EOF_ERA20C.mat', ...
    '~/Documents/BODELE/Modes/DATA/EOF_MERRA.mat', ...
    '~/Documents/BODELE/Modes/DATA/EOF_NCEP2.mat', ...
    '~/Documents/BODELE/Modes/DATA/EOF_NNRP.mat'};
names = {'CIRES-20CR','ERA-20CR','MERRA','NCEP2','NNRP'};


%% Plot SVD output
f1 = figure(1);
f1.Units = 'centimeters';
f1.PaperUnits = 'centimeters';
f1.Position = [0 0 18.3 20.7];
f1.PaperPosition = [0 2 18.3 20.7];
f1.PaperSize = [18.3 24.7];

cm = colormap('parula');
cm(1,:) = ones(1,3);
colormap(cm);

% let = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p', ...
%     'q','r','s','t','u','v','w','x','y','z'};
    
for i = 1:5
    
    load(files{i});
    [x_,y_] = meshgrid(lon,lat);
    x_ = x_'; y_ = y_';
    s = size(x_);
    x_ = reshape(x_,s(1),s(2));
    y_ = reshape(y_,s(1),s(2));
    
    ru_ = ru; rv_ = rv;
    
    % EOF Map
    for j = 1:3
        subplot(10,3,j+(i-1)*6)
        z = sqrt(ru(:,:,j).^2+rv(:,:,j).^2);
        ru(:,:,j) = ru(:,:,j)./z;
        rv(:,:,j) = rv(:,:,j)./z;
        z(isnan(z)) = 0;
        imagesc(lon,lat,z')
        cb = colorbar;
        set(cb,'fontsize',6)
        caxis([0 floor(max(z(:)))]);
        box on; hold on
        plot(coast.long,coast.lat,'k','linewidth',1);
        colormap(cm);

        di = 5;
        if lon(2)-lon(1) > .7; di = 3; end
        if lon(2)-lon(1) > 1; di = 2; end
        q = quiver(x_(1:di:end,1:di:end),y_(1:di:end,1:di:end), ...
            ru(1:di:end,1:di:end,j),rv(1:di:end,1:di:end,j),.5,'k');
        hold off
        
        % Text and position
        set(gca,'xlim',[-20 30],'ylim',[5 35]);
        set(gca,'ydir','normal','fontsize',6,'XAxisLocation','top');  
        if i==1;
            xlabel('Lon.')
            ttl = title(['EOF/PC ' num2str(j)]);
            ttl.FontWeight = 'bold';
            ttl.FontSize = 6;
        end
        if j==3; ylabel(cb,'Wspd.'); end
        if j==1; ylabel('Lat.'); end
        pos = get(gca,'position');
        pos(3) = 0.1140;
        set(gca,'position',pos+[0 0 0.04 0])

        vtxt = text(28,30,[num2str(round(100*sig(j)/sum(sig))) '%']);
        vtxt.BackgroundColor = [1 1 1];
        vtxt.EdgeColor = [0 0 0];
        vtxt.HorizontalAlignment = 'right';
        if j==1;
            txt = text(-40,5,names(i),'fontsize',6,'rotation',90);
            txt.HorizontalAlignment = 'center';
            txt.FontWeight = 'bold';
        end
        cb.Position = cb.Position-[0.01 0 0 0];
        caxis([0 floor(8*max(z(:)))/10 ])
    end
    
    for j = 1:3;
        subplot(10,3,j+3+(i-1)*6)
        pc = smooth(PC(:,j)./std(PC(:,j)),13);
        pc = (pc-mean(pc))/std(pc);
        plot(time(7:end-6),pc(7:end-6)) % uv
        hold on; 
        plot(time,time*0,'k')
        grid on; box on; hold off
        set(gca,'fontsize',6,'xlim',[min(time(7:end)) max(time(1:end-6))]);
        set(gca,'position',get(gca,'position')-[0 -.032 .015 .015])
        set(gca,'ylim',[floor(min(pc(7:end-6))) ceil(max(pc(7:end-6)))])
        if i==5; xlabel('Years'); end
        if j==1; ylabel('Std. dev.'); end
        if i==1; set(gca,'Xlim',[1850 max(time)]); end
        if i==1; set(gca,'Xlim',[1900 max(time)]); end
    end
    
    % write csv files
    [x,y] = meshgrid(lon,lat);
    lon = x'; lat = y';
    s = size(ru);
    A = [reshape(lon,s(1)*s(2),1) reshape(lat,s(1)*s(2),1) ...
        reshape(ru_(:,:,1),s(1)*s(2),1) reshape(rv_(:,:,1),s(1)*s(2),1) ...
        reshape(ru_(:,:,2),s(1)*s(2),1) reshape(rv_(:,:,2),s(1)*s(2),1) ...
        reshape(ru_(:,:,3),s(1)*s(2),1) reshape(rv_(:,:,3),s(1)*s(2),1) ];
    csvwrite(['figureE7.abc.' files{i}(35:end-4) '.csv'],A)
    
    st = size(time); if st(1)==1; time = time'; end
    A = [time PC(:,1:3)];
    csvwrite(['figureE7.def.' files{i}(35:end-4) '.csv'],A)

end


%% Print
print(f1,'-dpdf','-cmyk','figureE7.pdf')



%test


% figure(1)
% plot(Gau, ISREF,'ro')
% hold on
% plot(Gau, ISRWA,'x','Markersize',10)
% 
% xlim([min(Gau) max(Gau)])
% ylim([0 0.02])
% set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
%     'Courier','FontWeight','bold')
% xlabel('Gaussian noise (FWHM)')
% xlabh3 = get(gca,'XLabel');
% set(xlabh3,'Position',get(xlabh3,'Position') + [0 0 0]) 
% ylabel('ISR')
% hold off
% 
% 
% 
% figure(2)
% plot( COR,'x')
% hold on
% plot(Q./mean(Q),'rx')
% plot(PC./mean(PC),'gx')

%xlim([min(Gau) max(Gau)])
%ylim([0 0.015])
% set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
%     'Courier','FontWeight','bold')
% xlabel('Gaussian noise (FWHM)')
% xlabh3 = get(gca,'XLabel');
% set(xlabh3,'Position',get(xlabh3,'Position') + [0 0 0]) 
% ylabel('PC')

% figure(3)
% g = subplot(1,2,1);
% plot(DSTREAM.DATA.phase,DSTREAM.DATA.raw(:,1),'x')
% hold on
% plot(DSTREAM.DATA.phase,DSTREAM.DATA.NMODEL.signal(:,1)+1,'rx')
% set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
%     'Courier','FontWeight','bold')
% xlabel('Phase')
% ylabel('Rel. Flux')
% ax=get(g,'Position');
% ax(3)=ax(3)+0.1; %may need something better
% ax(1) = ax(1)-0.06;
% xlim([min(DSTREAM.DATA.phase) max(DSTREAM.DATA.phase)])
% ylim([0.992 1.008])
% set(g,'Position',ax, 'FontSize',15,'FontName','Courier','FontWeight','bold'); 
% hold off
% 
% 
% g = subplot(1,2,2);
% plot(DSTREAM.DATA.phase,DSTREAM.DATA.NMODEL.noise(:,1)+1,'x')
% 
% % set(gca,'ytick',[]);
% ax=get(g,'Position');
% ax(3)=ax(3)+0.1; %may need something better
% ax(1) = ax(1)-0.06;
% %xlim(min(waxvals), max(waxvals))
% %ylim(min(wayvals),max(wayvals))
% xlim([min(DSTREAM.DATA.phase) max(DSTREAM.DATA.phase)])
%  ylim([0.992 1.008])
% set(g,'Position',ax, 'FontSize',15,'FontName','Courier','FontWeight','bold')
% hold off


figure(4)
plot(DSTREAM.DATA.phase,DSTREAM.DATA.raw(:,1),'o','MarkerSize',5)
hold on
plot(DSTREAM.DATA.phase,DSTREAM.DATA.NMODEL.signal(:,1)+1,'rx')
plot(DSTREAM.DATA.phase,DSTREAM.DATA.NMODEL.noise(:,1)+1-0.008,'ks','MarkerSize',4)
hold off
xlim([min(DSTREAM.DATA.phase) max(DSTREAM.DATA.phase)])
ylim([0.99, 1.01])
set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
    'Courier','FontWeight','bold')
xlabel('Phase')
xlabh3 = get(gca,'XLabel');
set(xlabh3,'Position',get(xlabh3,'Position') + [0 0 0]) 
ylabel('Rel. flux')
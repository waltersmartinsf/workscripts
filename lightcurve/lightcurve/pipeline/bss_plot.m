%module plotting some of the bss_pipeline output

function bssp = bss_plot(DSTREAM,varargin)


% PHASE = DSTREAM.DATA.phase;
% minph = min(PHASE);
% maxph = max(PHASE);
% 
% 
% %plotting premix data if given
% if nargin == 2
%     PRE = varargin{1};
%     [ps1,ps2] = size(PRE);
%     bssp = figure(10);
% 
%     for i=1:ps2 
%         h = subplot(ps2,1,i);
%         plot(PHASE(:,1),PRE(:,i),'x','markersize',5)
%         ax=get(h,'Position');
%         ax(4)=ax(4)+0.042; %may need something better
%         set(h,'Position',ax); 
%         xlim([minph,maxph]);
%         set(gca,'xtick',[],'ytick',[]);
%     end
%     set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
%         'Courier','FontWeight','bold')
%     xlabel('Phase')
%     xlabh1 = get(gca,'XLabel');
%     set(xlabh1,'Position',get(xlabh1,'Position') + [0 0.01 0]) 
%     ylabel('Rel. flux')
% end
% 
% MIX = DSTREAM.DATA.raw;
% DEMIX = DSTREAM.DATA.MCOMBI.signals;
% 
% 
% %plotting mixed signals
% 
% [ms1,ms2] = size(MIX);
% figure(11)
% for i=1:ms2 
%         h = subplot(ms2,1,i);
%         plot(PHASE(:,1),MIX(:,i),'x','markersize',5)
%         ax=get(h,'Position');
%         ax(4)=ax(4)+0.042; %may need something better
%         set(h,'Position',ax); 
%         xlim([minph,maxph]);
%         set(gca,'xtick',[],'ytick',[]);
% end
% set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
%     'Courier','FontWeight','bold')
% xlabel('Phase')
% ylabel('Rel. flux')
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [0.0015 0 0]) 
% xlabh2 = get(gca,'XLabel');
% set(xlabh2,'Position',get(xlabh2,'Position') + [0 0.001 0]) 
% 
% %plotting demixed signals unsorted (in case figure13 fucks up)
% 
% figure(12)
% for i=1:ms2 
%         h = subplot(ms2,1,i);
%         plot(PHASE(:,1),DEMIX(i,:),'x','markersize',5)
%         ax=get(h,'Position');
%         ax(4)=ax(4)+0.042; %may need something better
%         set(h,'Position',ax); 
%         xlim([minph,maxph]);
%         set(gca,'xtick',[],'ytick',[]);
% end
% set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
%     'Courier','FontWeight','bold')
% xlabel('Phase')
% xlabh3 = get(gca,'XLabel');
% set(xlabh3,'Position',get(xlabh3,'Position') + [0 0.8 0]) 
% ylabel('Rel. flux')
% 
% %plotting demixed signals sorted after noise, systematic noise is red
% 
% [desn1,desn2] = size(DSTREAM.DATA.FIND.sn);
% [dewn1,dewn2] = size(DSTREAM.DATA.FIND.wn);
% [delc1,delc2] = size(DSTREAM.DATA.FIND.lc);
% figure(13)
% for i=1:desn2 
%         h = subplot(ms2,1,i);
%         plot(PHASE(:,1),DSTREAM.DATA.FIND.sn(:,i),'rx','markersize',5)
%         ax=get(h,'Position');
%         ax(4)=ax(4)+0.042; %may need something better
%         set(h,'Position',ax); 
%         set(gca,'xtick',[],'ytick',[]);
%         xlim([minph,maxph]);
% end
% for i=1:dewn2 
%         h = subplot(ms2,1,i+desn2);
%         plot(PHASE(:,1),DSTREAM.DATA.FIND.wn(:,i),'x','markersize',5)
%         ax=get(h,'Position');
%         ax(4)=ax(4)+0.042; %may need something better
%         set(h,'Position',ax); 
%         set(gca,'xtick',[],'ytick',[]);
%         xlim([minph,maxph]);
% end
% for i=1:delc2 
%         h = subplot(ms2,1,i+desn2+dewn2);
%         plot(PHASE(:,1),DSTREAM.DATA.FIND.lc(:,i),'x','markersize',5)
%         ax=get(h,'Position');
%         ax(4)=ax(4)+0.042; %may need something better
%         set(h,'Position',ax); 
%         set(gca,'xtick',[],'ytick',[]);
%         xlim([minph,maxph]);
% end
% 
% set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
%     'Courier','FontWeight','bold')
% xlabel('Phase')
% xlabh = get(gca,'XLabel');
% set(xlabh,'Position',get(xlabh,'Position') + [0 0.8 0]) 
% ylabel('Rel. flux')


%plotting hinton diagrams for ISR_ef and ISR_ef matrices 

w = DSTREAM.INFO.MCOMBI.ISRef1;
w2 = DSTREAM.INFO.MCOMBI.ISRwa1;
scale1 = max(max(abs(w)));
scale2 = max(max(abs(w2)));
diffscale = scale2 ./scale1

if diffscale <= 1.0
    [efxvals, efyvals, efcolor] = hintmat(DSTREAM.INFO.MCOMBI.ISRef1,1);
    [waxvals, wayvals, wacolor] = hintmat(DSTREAM.INFO.MCOMBI.ISRwa1,diffscale);
elseif diffscale > 1.0
     diffscale = scale1 ./ scale2;
     [efxvals, efyvals, efcolor] = hintmat(DSTREAM.INFO.MCOMBI.ISRef1,diffscale);
     [waxvals, wayvals, wacolor] = hintmat(DSTREAM.INFO.MCOMBI.ISRwa1,1);
end






figure(15)
set(gcf, 'Units', 'pixels', 'Position', [0 0 700 400])
g = subplot(1,2,1);
patch(efxvals', efyvals', efcolor', 'Edgecolor', 'none');
%set(gca,'xtick',[1,2,3,4,5],'ytick',[1,2,3,4,5]);
ax=get(g,'Position');
ax(3)=ax(3)+0.1; %may need something better
ax(1) = ax(1)-0.06;
%xlim(min(efxvals),max(efxvals))
%ylim(min(efyvals),max(efyvals))
set(g,'Position',ax, 'FontSize',15,'FontName','Courier','FontWeight','bold'); 
title('ISR-EFICA', 'FontSize',15,'FontName','Courier','FontWeight','bold')
text(0.5,0.04,['max: ' num2str(max(max(abs(DSTREAM.INFO.MCOMBI.ISRef1))),5)],...
    'Units','normalized','FontSize',15,'FontName','Courier','FontWeight','bold')


g=subplot(1,2,2);
patch(waxvals', wayvals', wacolor', 'Edgecolor', 'none');
%set(gca,'xtick',[1,2,3,4,5],'ytick',[]);
set(gca,'ytick',[]);
ax=get(g,'Position');
ax(3)=ax(3)+0.1; %may need something better
ax(1) = ax(1)-0.06;
%xlim(min(waxvals), max(waxvals))
%ylim(min(wayvals),max(wayvals))
set(g,'Position',ax, 'FontSize',15,'FontName','Courier','FontWeight','bold'); 
title('ISR-WASOBI', 'FontSize',15,'FontName','Courier','FontWeight','bold')
text(0.5,0.04,['max: ' num2str(max(max(abs(DSTREAM.INFO.MCOMBI.ISRwa1))),5)],...
    'Units','normalized','FontSize',15,'FontName','Courier','FontWeight','bold')
hold off



%plotting fitted data if existing 

















%bss post processing (so far only autcorrelation check)

function [DSTREAM] = bss_post(DSTREAM)

% resid = DSTREAM.DATA.FIT.resid(120:240);
resid = DSTREAM.DATA.raw(120:240,8);


resid = resid(:,1);

[ACF,lags,bounds] = autocorr(resid,30,0,2);

figure17 = figure;
% Create axes
axes1 = axes('Parent',figure17);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 100]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-0.8 1]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

% Create stem
stem(lags,ACF,'MarkerSize',4,'MarkerFaceColor','auto','Color',[1 0 0],'linewidth',1);
hold on
hline(bounds(1),'blue')
hline(bounds(2),'blue')
set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName','Courier','FontWeight','bold')
xlabel('Lags')
ylim([-1 1])
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [0 0 0]) 
ylabel('Autocorrelation')
%bss pipeline quick lightcurve plotting routine 


function DSTREAM = bss_quickplot(DSTREAM, LCNO)


figure(10)
plot(DSTREAM.DATA.raw(:,LCNO),'s', 'color',[0.6 0.6 0.6],'MarkerSize',4)
hold on
plot(DSTREAM.DATA.NMODEL.signal(:,LCNO),'x')
plot(DSTREAM.DATA.NMODEL.noise(:,LCNO)+1,'rx')
hold off

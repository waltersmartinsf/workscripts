%second version of the noise-model/lightcurve fitting algorithm



function [DSTREAM] = bss_fit2(DSTREAM,PSTREAM,varargin)

% global DSTREAM
% global PSTREAM

if nargin == 3
    DSTREAM.INFO.FIT.lcindex = varargin{1};
else
    DSTREAM.INFO.FIT.lcindex = PSTREAM.FIT.lcindex;
end

%getting raw data vectors
RAW = DSTREAM.DATA.raw(:,DSTREAM.INFO.FIT.lcindex);
RAWmean = mean(RAW);
RAW = RAW - RAWmean;
[RAWs1, RAWs2] = size(RAW);
if RAWs1 < RAWs2
    RAW = transpose(RAW);
end

%getting systematic noise components 
SN = DSTREAM.DATA.FIND.sn;
% SN = cat(2,SN,(PRE.lc-mean(PRE.lc)));

[SNs1,SNs2] = size(SN);

for i=1:SNs2
    SN(:,i) = SN(:,i) - mean(SN(:,i));
end

PHASE = DSTREAM.DATA.phase;
% minph = min(PHASE);
% maxph = max(PHASE);

%constructing SN eigenvalue vector using unity priors
lcpar0 = ones(SNs2+1,1);

%starting minimization over SNeig
options.MaxFunEvals = 100000;
[lcpar1] = fminsearch(@(lcpar0) lcss2(lcpar0,PHASE,RAW,SN),lcpar0,options);


%constructing models, fits, etc 
mix = zeros(size(SN));

for i=1:SNs2;
    mix(:,i) = lcpar1(i,1) .* SN(:,i);
end
DSTREAM.DATA.FIT.vectors = mix;

mixadd = sum(mix,2);
resid = RAW - mixadd;
lcpar1 = diag(lcpar1);

%loading results into DSTREAM
DSTREAM.DATA.NMODEL.nscale = lcpar1;
DSTREAM.DATA.FIT.model = mixadd;
DSTREAM.DATA.FIT.resid = resid;
%DSTREAM.DATA.FIT.agol = mamod;


% %plotting result (will be moved into bss_plot soon)
% figure(20)
% clf
% plot(PHASE, RAW+1, 'x')
% hold on
% plot(PHASE,mixadd+1, 'r','linewidth',1.0)
% plot(PHASE,resid+0.98,'kx')
% xlim([minph,maxph]);
% %ylim([0.975, 1.01])
% set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
%     'Courier','FontWeight','bold')
% xlabel('Phase')
% xlabh3 = get(gca,'XLabel');
% set(xlabh3,'Position',get(xlabh3,'Position') + [0 0 0]) 
% ylabel('Rel. flux')
% hold off


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ss] = lcss2(SNeig,PHASE,RAW,SN)
global DSTREAM
global PSTREAM


[SNs1,SNs2] = size(SN);

%setting up mixing array
mix = zeros(SNs1,SNs2+1);
for i=1:SNs2;
    mix(:,i) = SNeig(i,1) .* SN(:,i);
end


%calculating Mandel&Agol lightcurve as last SNeig entry
rs = PSTREAM.FIT.LCPARAMS.rstar;
inc = PSTREAM.FIT.LCPARAMS.inc;
sma = PSTREAM.FIT.LCPARAMS.sma;
u1 = PSTREAM.FIT.LCPARAMS.u1;
u2 = PSTREAM.FIT.LCPARAMS.u2;
ecc = PSTREAM.FIT.LCPARAMS.ecc;
omega = PSTREAM.FIT.LCPARAMS.omega;

[SNeigs1,SNeigs2] = size(SNeig);

p = SNeig(SNeigs2);

mamod = agolquad2(p,rs,inc,sma,PHASE,u1,u2,ecc,omega);
mamod = mamod - mean(mamod);


%adding model to mix now 
mix(:,SNs2+1) = mamod;


mixadd = sum(mix,2);

ss = sum(sum((RAW -(mixadd)).^2));

figure(1)
plot(RAW+1,'x')
hold on
plot(mixadd+1,'r','linewidth',1.0)
hold off

end



function [FRlimb,FRnolimb] = agolquad2(p,radstar,inclination,semimajoraxis,phase,u1,u2,ecc,omega)


% radplanet = 1.184;
%p = 0.1473;
p = 0.126;
    
%converting units to cgs
inclination = inclination / 360 * 2*pi;
semimajoraxis = semimajoraxis *149.6e10;
radstar = radstar * 6.955e9;
%  radplanet = radplanet * 7.1492e8;


%  rp   = double(radplanet);
% rs   = double(radstar);
%  p    = double(rp/rs);
u1   = double(u1);
u2   = double(u2);

% go from phase to center-to-center distance
%distancevector = delta((phase-0.5d0)*2d0*!dpi,inclination)*semimajoraxis/radstar

distancevector = delta(((phase)*2d0*pi),inclination,ecc,omega)*semimajoraxis/radstar;
z    = double(distancevector);

% executing Mandel & Agol 2002
[FRlimb,FRnolimb] = agolquad_sub(z,u1,u2,p);

end
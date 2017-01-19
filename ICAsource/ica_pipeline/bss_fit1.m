%first version of the noise-model/lightcurve fitting algorithm


function [DSTREAM] = bss_fit1(LC,varargin)

global DSTREAM;
global PSTREAM;

if nargin > 1
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
SN = cat(2,SN,(LC-mean(LC)));


[SNs1,SNs2] = size(SN);

for i=1:SNs2
    SN(:,i) = SN(:,i) - mean(SN(:,i));
end

%LC = PREMIX.lc;
PHASE = DSTREAM.DATA.phase;
minph = min(PHASE);
maxph = max(PHASE);

%doing white-noise reduction 
% ksr = bss_ksr(PHASE,RAW(:,1),0.0001);
% raw1 = ksr.f;


%constructing SN eigenvalue vector using unity priors
[SNs1,SNs2] = size(SN);
SNeig = ones(SNs2,1);
%SNeig = diag(SNeig)

size(SNeig);

%starting minimization over SNeig


fig1 = figure(1);
% global winsize 
% winsize = get(fig1,'Position');
% numframes = 1300;
% global A
% A = moviein(numframes,fig1,winsize);

plot(RAW,'x')
set(gca,'FontSize',20)
xlabel('No. of spectra','FontSize',20)
ylabel('Relative flux','FontSize',20)
% 
% global ind
% ind = 1;
% A(ind) = getframe(fig1);


options.MaxFunEvals = 1000000;

[lcpar0] = fminsearch(@lcss,SNeig,options);

mix = zeros(size(SN));
for i=1:SNs2;
    mix(:,i) = lcpar0(i,1) .* SN(:,i);
end
mixadd = sum(mix,2);
resid = RAW - mixadd;
lcpar0 = diag(lcpar0);

%loading results into DSTREAM
DSTREAM.DATA.NMODEL.nscale = lcpar0;
DSTREAM.DATA.FIT.model = mixadd;
DSTREAM.DATA.FIT.resid = resid;


%plotting result (will be moved into bss_plot soon)
figure(20)
clf
plot(PHASE, RAW+1, 'x')
hold on
plot(PHASE,mixadd+1, 'ro','linewidth',1)
plot(PHASE,resid+0.98,'kx')
xlim([minph,maxph]);
ylim([0.975, 1.015])
set(gca,'xtickMode', 'auto','ytickMode','auto', 'FontSize',15,'FontName',...
    'Courier','FontWeight','bold')
xlabel('Phase')
xlabh3 = get(gca,'XLabel');
set(xlabh3,'Position',get(xlabh3,'Position') + [0 0 0]) 
ylabel('Rel. flux')

end


function [ss] = lcss(SNeig)

global DSTREAM
global PRE

%getting systematic noise components 
SN = DSTREAM.DATA.FIND.sn;
%SN = cat(2,SN,DSTREAM.DATA.FIND.lc);
SN = cat(2,SN,(PRE.lc-mean(PRE.lc)));
[SNs1,SNs2] = size(SN);

% %getting raw data vectors
% RAW = DSTREAM.DATA.raw(:,1);
% RAW = RAW - mean(RAW);
% [RAWs1, RAWs2] = size(RAW);
% if RAWs1 < RAWs2
%     RAW = transpose(RAW);
% end


%getting raw data vectors
RAW = DSTREAM.DATA.raw(:,DSTREAM.INFO.FIT.lcindex);
RAWmean = mean(RAW);
RAW = RAW - RAWmean;
[RAWs1, RAWs2] = size(RAW);
if RAWs1 < RAWs2
    RAW = transpose(RAW);
end


mix = zeros(size(SN));

size(SNeig);
[SNs1,SNs2] = size(SN);

for i=1:SNs2
    SN(:,i) = SN(:,i) - mean(SN(:,i));
end



for i=1:SNs2;
    mix(:,i) = SNeig(i,1) .* SN(:,i);
end


%mix = SNeig * SN';
mixadd = sum(mix,2);

size(mix);
size(mixadd);
size(RAW);

ss = sum((RAW -(mixadd)).^2);
ss = sum(ss);


% figure(3)
% plot(RAW' - mixadd,'x')

fig1 = figure(1);
% winsize = get(fig1,'Pplotosition');
plot(RAW,'x')
hold on
plot(mixadd,'r','linewidth',1)
hold off
set(gca,'FontSize',20)
xlabel('No. of spectra','FontSize',20)
ylabel('Relative flux','FontSize',20)

% global ind
% ind = ind +1
% global A
% A(ind) = getframe(fig1);
% 
% global BLE;
% BLE.A = A; 
end

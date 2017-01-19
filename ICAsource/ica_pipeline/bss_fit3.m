%bss_fit3 - optimised for HST/NICMOS HD189 primary eclipse data 
%it assumes a full fledged DSTREAM/PSTREAM structure as input 
%it does NOT model a Mandel & Agol lightcurve like bss_fit2 but only
%computes the noisemodel fit on the out-of-transit orbits. 


function [DSTREAM] = bss_fit3(DSTREAM)

OOTind1 = 128;      %345;%247;
OOTind2 = 262;      %476;%378;

ITexclude = OOTind1:OOTind2;
DSTREAM.INFO.FIT.intranind = ITexclude;


%getting raw data vectors
RAW = DSTREAM.DATA.raw;
[RAWs1,RAWs2] = size(RAW);

OOTnewsize = RAWs1 - (OOTind2-OOTind1)-1;

%getting systematic noise components 
SN = DSTREAM.DATA.FIND.sn;
[SNs1,SNs2] = size(SN);

%trimming systematics 
SNcut = zeros(OOTnewsize,SNs2);
for i=1:SNs2
    cuttmp = SN(:,i);
    cuttmp(ITexclude) = [];
    SNcut(:,i) = cuttmp;
end

%normalising systematics to zero mean
for i=1:SNs2
    SNcut(:,i) = SNcut(:,i) - mean(SNcut(:,i));
end

lcparfin = zeros(SNs2,RAWs2);


%loop over all RAWs2 channels
for i=1:RAWs2
    
    data = RAW(:,i);
    
    %trimming data
    datacut = data;
    datacut(ITexclude) =[];
    
    %normalising data to zero mean
    datacutmean = mean(datacut);
    datacut = datacut - datacutmean;
    
    %constructing SN eigenvalue vector using unity priors
    lcpar0 = ones(SNs2,1);
    
    %starting minimization over SNeig
    options.MaxFunEvals = 100000;
    [lcpar1] = fminsearch(@(lcpar0) lcss3(lcpar0,datacut,SNcut),lcpar0,options);
    
    lcparfin(:,i) = lcpar1;
    
end


%normalising systematics to zero mean
for i=1:SNs2
    SN(:,i) = SN(:,i) - mean(SN(:,i));
end


%constructing models, fits, etc 
mixtotal = zeros(RAWs1,RAWs2);

for i=1:RAWs2
    mix = zeros(SNs1,SNs2);
    for j=1:SNs2;
        mix(:,j) = lcparfin(j,i) .* SN(:,j);
    end
    mixadd = sum(mix,2);

    mixtotal(:,i) = mixadd;
end

%displaying every fit
% figure(1)
% plot(datacut,'x')
% hold on
% plot(mixadd,'r','linewidth',1.2)
% hold off


% DSTREAM.DATA.FIT.vectors = mix;

resid = RAW - mixtotal;
residmean = mean(resid,2);
residmeanoot = residmean; 
residmeanoot(ITexclude) = [];
off = mean(residmeanoot) - 1;

residmean = residmean - off;
resid = resid - off;


% lcparfin = diag(lcparfin);

%loading results into DSTREAM
DSTREAM.DATA.NMODEL.nscale = lcparfin;
DSTREAM.DATA.FIT.model = mixtotal + 1;
DSTREAM.DATA.FIT.resid = resid;
% DSTREAM.DATA.FIT.trandep = 1-inmean;
%DSTREAM.DATA.FIT.agol = mamod;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ss] = lcss3(SNeig,DATA,SN)

[SNs1,SNs2] = size(SN);

%setting up mixing array
mix = zeros(SNs1,SNs2);
for i=1:SNs2;
    mix(:,i) = SNeig(i,1) .* SN(:,i);
end


mixadd = sum(mix,2);

ss = sum(sum((DATA -(mixadd)).^2));

% figure(1)
% plot(DATA,'x')
% hold on
% plot(mixadd,'r','linewidth',1.2)
% hold off

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%agolquad is maintained here but not used.

function [FRlimb,FRnolimb] = agolquad3(p,radstar,inclination,semimajoraxis,phase,u1,u2,ecc,omega)


inclination = 85.51;
semimajoraxis = 0.03142;
radstar = 0.788;
u1 = 0.0104;
u2 = 0.4663;
    


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
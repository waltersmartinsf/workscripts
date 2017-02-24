%module reconstructing the appropriate noise-model for each input channel
%it requires inputs from bss_find and bss_mcombi

function [DSTREAM,PSTREAM] = bss_nmodel(DSTREAM,PSTREAM)


W = DSTREAM.DATA.MCOMBI.w; 
S = DSTREAM.DATA.MCOMBI.signals;
sn = DSTREAM.DATA.FIND.INDEX.sn;
lc = DSTREAM.DATA.FIND.INDEX.lc;

[Ws1,Ws2] = size(W);

%recreating the full mixing matrix where W = A^(-1)
A = inv(W);

%creating systematic noise only model
NOISE = A(:,sn) * S(sn,:);
NOISE = transpose(NOISE);

%creating pure signal only model
SIGNAL = A(:,lc) * S(lc,:);
SIGNAL = transpose(SIGNAL);

%adding original mean to SIGNAL
[Ss1,Ss2] = size(SIGNAL);

for i=1:Ss2
    SIGNAL(:,i) = SIGNAL(:,i) + DSTREAM.DATA.mean(i);
end


%creating scaling coefficient vector
NMODSCALE = ones(1,Ws2);
DIAG = diag(NMODSCALE);

%adding data to data stream
DSTREAM.DATA.NMODEL.noise = NOISE;
DSTREAM.DATA.NMODEL.signal = SIGNAL;
DSTREAM.DATA.NMODEL.a = A;
DSTREAM.DATA.NMODEL.nscale = DIAG;
DSTREAM.INFO.NMODEL.nscalefit = 0;









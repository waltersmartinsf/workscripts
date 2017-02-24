%module performing a gaussian kernel regression smoothing for white noise
%reduction 


function [DSTREAM, PSTREAM] = bss_ksr(DSTREAM,PSTREAM)

bandwidth = PSTREAM.KERNEL.bandwidth;
PHASE = DSTREAM.DATA.phase;
DATA = DSTREAM.DATA.raw;
[Ds1, Ds2] = size(DATA);

DATA2 = zeros(size(DATA));

for i=1:Ds2
    r = ksr(PHASE(:,1),DATA(:,i),bandwidth,length(DATA));
    DATA2(:,i) = real(r.f);
end

DSTREAM.DATA.raw_orig = DSTREAM.DATA.raw;
DSTREAM.DATA.raw = DATA2;

%  figure(5)
%  plot(PHASE(:,1),DATA(:,1),'x',r.x,r.f,'r--','linewidth',2)


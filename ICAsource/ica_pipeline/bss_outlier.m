% outlier rejection


function [OUT,PHASEOUT] = bss_outlier(DATA,PHASE,BLOCK,sigma)

[s1,s2] = size(DATA);
if s1 < s2
    DATA = DATA';
end
[s1,s2] = size(DATA);


%sigma = 1;

L = 1;
H = BLOCK;

OUT = zeros(0,1);
PHASEOUT = zeros(0,1);



for i=1:s1
    tmp = DATA(L:H,1);
    phasetmp = PHASE(L:H,1);
    tmpmean = mean(tmp);
    tmpstd = std(tmp);
    
    for i=1:length(tmp)
        if (tmp(i) < tmpmean+tmpstd*sigma) && (tmp(i) > tmpmean-tmpstd*sigma)
            OUT = cat(1,OUT,tmp(i));
            PHASEOUT = cat(1,PHASEOUT,phasetmp(i));
        end
    end
    H = H + BLOCK;
    L = L + BLOCK;
    if (H >= s1),break;end
end



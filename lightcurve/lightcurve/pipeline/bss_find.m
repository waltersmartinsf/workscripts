%signal identification algortihm. It compared components obtained from
%bss_mcombi to the first PC obtained from bss_pca to identify the
%lightcurve signal. Other signals will be judged by a ljung-box test and if
%significantly non-guassian, accepted for the systematic noise model. 
%
%To be added:
%
%-if too many white noise components found, change PSTREAM value to force a
%PCA dimensionality reduction on the next run 
%-identify what to do when more than one LC is found. Build in a manual
%inspection and rejection mechanism
%-many things that I cannot think of right now


function [DSTREAM,PSTREAM] = bss_find(DSTREAM,PSTREAM)

%loading principal and mcombi components from data-stream
PC = DSTREAM.DATA.PCA.score;
FiPC = PC(:,1);
%[PCs1, PCs2] = size(PC);

IC = DSTREAM.DATA.MCOMBI.signals';
[ICs1,ICs2] = size(IC);

%checking data stream for any nonseparated componenents
if exist('DSTREAM.MCOMBI.nonsep','var') ~= 0 
    fprintf('WARNING: non-separated components detected! \n')
    DSTREAM.INFO.WARNINGS.find1 = 'WARNING: non-separated components detected!';
end


%separating IC containing lightcurve signal from the others 
INDEX = zeros(1,ICs2);
for i=1:ICs2;
    INDEX(1,i) = i;
end
INDEX_lc = zeros(1,0);
INDEX_nonlc = zeros(1,0);
INDEX_sn = zeros(1,0);
INDEX_wn = zeros(1,0);

nLC = 0;
LC = zeros(length(IC),0);
nonLC = zeros(length(IC),0);
corrstats = zeros(2,0);
for i=1:ICs2
    [rho, pval] = corr(IC(:,i),FiPC);
    tmp = zeros(2,1);
    tmp(1,1) = rho;
    tmp(1,2) = pval;
    corrstats = cat(2,corrstats,tmp);
    rho2 = sqrt(rho^2);
    if rho2 > PSTREAM.FIND.pearsonthr
        nLC = nLC + 1;
        LC = cat(2,LC,IC(:,i));
        INDEX_lc = cat(2,INDEX_lc,INDEX(1,i));
    else
        nonLC = cat(2,nonLC,IC(:,i));
        INDEX_nonlc = cat(2,INDEX_nonlc,INDEX(1,i));
        
    end
end
if nLC > 1
    fprintf('WARNING: more than one light-curve signal component detected! \n')
    DSTREAM.INFO.WARNINGS.find2 = 'WARNING: more than one light-curve signal component detected! \n';
end


%clasifying remainign ICs according to their white or non-white
%autocovariance matrices using Ljung-Box test

[nonLCs1, nonLCs2] = size(nonLC);

WN = zeros(length(nonLC),0);
SN = zeros(length(nonLC),0);
LBQstats = zeros(4,0);
WNcount = 0;
SNcount = 0;
for i=1:nonLCs2
    [h,pValue,stat,cValue] = lbqtest(nonLC(:,i),'alpha',PSTREAM.LBQ.alpha);
    tmp = zeros(4,1);
    tmp(1,1) = h;
    tmp(2,1) = pValue;
    tmp(3,1) = stat;
    tmp(4,1) = cValue;
    LBQstats = cat(2,LBQstats,tmp);
   
    if h == 1
        SN = cat(2,SN,nonLC(:,i));
        SNcount = SNcount +1;
        INDEX_sn = cat(2,INDEX_sn,INDEX_nonlc(1,i));
    elseif h == 0
        WN = cat(2,WN,nonLC(:,i));
        WNcount = WNcount + 1;
        INDEX_wn = cat(2,INDEX_wn,INDEX_nonlc(1,i));
    end
end
if WNcount > 1 
    fprintf('WARNING: number of Gaussian components in data: ')
    WNcount 
    fprintf('\n')
    DSTREAM.INFO.WARNINGS.find3 = 'WARNING: more than one Gaussian component in data';
end

%adding results to DSTREAM
DSTREAM.INFO.FIND.corrstats = corrstats;
DSTREAM.INFO.FIND.LBQstats = LBQstats;
DSTREAM.DATA.FIND.lc = LC;
DSTREAM.DATA.FIND.nonlc = nonLC;
DSTREAM.DATA.FIND.sn = SN;
DSTREAM.DATA.FIND.wn = WN;
DSTREAM.DATA.FIND.INDEX.lc = INDEX_lc;
DSTREAM.DATA.FIND.INDEX.nonlc = INDEX_nonlc;
DSTREAM.DATA.FIND.INDEX.sn = INDEX_sn;
DSTREAM.DATA.FIND.INDEX.wn = INDEX_wn;


    



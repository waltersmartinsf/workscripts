%function returning sum of squares for bss_fit1


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

figure(1)
plot(RAW,'x')
hold on
plot(mixadd,'r','linewidth',1)
hold off

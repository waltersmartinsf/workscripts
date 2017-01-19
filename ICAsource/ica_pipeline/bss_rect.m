% small code to rectify lightcurves with second order polynomial fits 
% to out of transit data


function [OUT] = bss_rect(LC,PHASE,LOW,UP)


[s1,s2] = size(LC);

if s2 > s1
    LC = transpose(LC);
end

%rectifying data
[lcs1,lcs2] = size(LC);

NOLC = zeros(0,1);
NOLC = cat(1,NOLC, LC(1:LOW,:));
NOLC = cat(1,NOLC,LC(UP:lcs1,:));

NOLCP = zeros(0,1);
NOLCP = cat(1,NOLCP, PHASE(1:LOW,:));
NOLCP = cat(1,NOLCP,PHASE(UP:lcs1,:));

p = polyfit(NOLCP(:,1),NOLC(:,1),1);
fit = polyval(p,PHASE(:,1));
fit2 = polyval(p,NOLC(:,1));

figure(15)
clf
plot(PHASE,LC,'x')
hold on 
plot(NOLCP,NOLC,'ro')
plot(PHASE,fit,'k')
hold off

OUT = LC - fit +1;
 
figure(16)
clf
plot(PHASE,OUT,'x')
hold off


% routine calculating upsilon test statistics from Waldmann & Tinetti 2011


function [Ups] = upsentr(D1,D2) 

[M1,N1] = size(D1)
[M2,N2] = size(D2);
%[M3,N3] = size(D3);

% [n1,xble] = hist(D1,100);
% bar(n1)
% [n2,xble] = hist(D2,100);
% %[n3,xble] = hist(D3,(M3./50));
% 
% Q1 = n1/M1;
% Q2 = n2/M2;
%Q3 = n3/M3;


Q1 = ksdensity(D1,'function','pdf');
Q2 = ksdensity(D2,'function','pdf');

Q1 = Q1 ./ sum(Q1);
Q2 = Q2 ./ sum(Q2);

figure(3)
plot(Q1)
hold on
plot(Q2,'r')



X = zeros(length(Q1),1);
for i=1:length(Q1);
    X(i) = (i-(length(Q1)/2))*1/length(Q1);
end

G = normpdf(X,0,0.12);
G = G ./ sum(G);

Q1 = Q1';
Q2 = Q2';

sizeq1 = size(Q1)
sizeq2 = size(Q2)
sizeg = size(G)
sumg = sum(G)
plot(G,'green')
hold off
legend('Q1','Q2','G')

c = 0.000000000001;

% DQ1 = - sum(G(:,1).*log2(G(:,1)+eps)) - sum(G(:,1).*log2(Q1(:,1)+eps))
% DQ2 = - sum(G(:,1).*log2(G(:,1)+eps)) - sum(G(:,1).*log2(Q2(:,1)+eps))
DQ1 = sum(Q1(:,1) .*(log2(Q1(:,1)+c)-log2(G(:,1)+c)))
DQ2 = sum(Q2(:,1) .*(log2(Q2(:,1)+c)-log2(G(:,1)+c)))
% DQ3 = sum(Q2(:,1) .*(log2(Q2(:,1))-log2(Q1(:,1))))

HG = -sum(G(:,1).*log2(G(:,1)+eps))
HQ1 = -sum(Q1(:,1).*log2(Q1(:,1)+eps))
HQ2 = -sum(Q2(:,1).*log2(Q2(:,1)+eps))


Ups = double(1. - (DQ2 ./ DQ1))


figure(4)
plot(Q1(:,1) .*(log2(Q1(:,1))-log2(G(:,1))))
hold on 
plot(Q2(:,1) .*(log2(Q2(:,1))-log2(G(:,1))),'red')
hold off




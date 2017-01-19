[% calculates shannon information entropy

function [H] = bss_shan(DATA)

[M,N] = size(DATA);

[n, xble] = hist(DATA,100);

norm = n/M;

H = -sum(norm(:,1).*log2(norm(:,1)+eps))

% figure(1)
% bar(norm)
%module doing PCA analysis, dimensionality reduction and PC signal
%extractino


function [DSTREAM, PSTREAM] = bss_pca(DSTREAM,PSTREAM)

%transposing data
[RAWs1,RAWs2] = size(DSTREAM.DATA.raw);
if RAWs1 < RAWs2
  RAWtrans = transpose(DSTREAM.DATA.raw);
else
    RAWtrans = DSTREAM.DATA.raw;
end

%measuring the mean and recording it
[Rs1,Rs2] = size(RAWtrans);
RAWMEAN = zeros(Rs2,1);

for i=1:Rs2
    RAWMEAN(i) = mean(RAWtrans(:,i));
end

DSTREAM.DATA.mean = RAWMEAN;

    
%%%%%%% IMPORTANT! TEMPORARILY DISABLED DUE TO UNEXPLAINED 'NAN' ERROR %%%

%doing PCA decomposition
[PC, SCORE, LATENT] = princomp(RAWtrans);

%saving outcome
DSTREAM.DATA.PCA.pc = PC;
DSTREAM.DATA.PCA.score = SCORE;
DSTREAM.DATA.PCA.latent = LATENT;
% 
% %performing dimensionality reduction if requested by PSTREAM. 
% if PSTREAM.PCA.maxcoeff > 0
%     PSTREAM.PCA.maxcoeff
%     DSTREAM.DATA.pca = SCORE(:,1:PSTREAM.PCA.maxcoeff) * PC(:,1:PSTREAM.PCA.maxcoeff)';
%     DSTREAM.INFO.PCA.pcretained = PSTREAM.PCA.maxcoeff;
% elseif PSTREAM.PCA.maxcoeff < 0
%     varretained = cumsum(LATENT)./sum(LATENT);
%     for i=1:length(varretained)
%         if varretained(i) >= PSTREAM.PCA.maxvar;
%             DSTREAM.INFO.PCA.pcretained = i;
%             DSTREAM.INFO.PCA.varretained = varretained(i);
%             break
%         end
%     end
%     DSTREAM.DATA.pca = SCORE(:,1:DSTREAM.INFO.PCA.pcretained) * PC(:,1:DSTREAM.INFO.PCA.pcretained)';
% elseif PSTREAM.PCA.maxcoeff == 0 
% %     DSTREAM.DATA.pca = SCORE * PC';
% %     DSTREAM.INFO.PCA.pcretained = PSTREAM.PCA.maxcoeff;
%    DSTREAM.DATA.pca = RAWtrans;

DSTREAM.DATA.pca = RAWtrans;

end


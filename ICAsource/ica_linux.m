% Here we will add all directories with matlab codes that are usefull 
%to reduce data
%addpath('/Users/walter/Dropbox/work/PythonUtility/ica_pipeline/')
addpath('/home/walter/GitHub/workscripts/ICAsource/ica_pipeline')

%let's import our data
%change for the new directory
%cd('/Users/walterwsmf/Dropbox/research/ica_project/xo2b/ica_tables/B_filter/')
%importing csv file remove the first row and first column (index values)
X = csvread('data_to_ica.csv',1,1)

%Obtain the whitening information
Xmean = mean(X')
Xstd  = std(X')

%whitening our data
%X = bsxfun(@minus, X, Xmean')
%X = bsxfun(@rdivide, X, Xstd')

%display our data
% plot(X','DisplayName','X')
% title('Mixing Signals')
% saveas(gca,'mixing_signals','png')

%Applied ICA-clustering in the data:
[result, ISRefica, ISRwasobi, Wmatrix, Wefica, Wwasobi, NoiseModel, nonseparablecomponents] = icapy(X)

%re-scale our signal
%signal = bsxfun(@times,Xstd',result)
%signal = bsxfun(@plus,Xmean',signal)

%display the results

%signals
% figure()
% plot(result','DisplayName','result')
% title('Reulsts signals from ICA')
% saveas(gca,'original_signals_result','png')

% figure()
% plot(signal','DisplayName','signal')
% title('Reulsts signals rescaled from ICA')
% saveas(gca,'original_signals_scaled','png')

Norm_ISR = (mean(ISRefica(:))+mean(ISRwasobi(:)))/2.

%display ISR results
% figure()
% image(ISRefica/Norm_ISR,'CDataMapping','scaled')
% set(gca,'YDir','normal')
% title(sprintf('ISR EFICA Normalized at %f',Norm_ISR))
% colorbar
% saveas(gca,'ISR_normalized_efica','png')

% %display ISR results
% figure()
% image(ISRwasobi/Norm_ISR,'CDataMapping','scaled')
% set(gca,'YDir','normal')
% title(sprintf('ISR WASOBI Normalized at %f',Norm_ISR))
% colorbar
% saveas(gca,'ISR_normalized_wasobi','png')

%comparing initial and final lightcurves


%Saving our results
csvwrite('ica_resuts.csv',result)

%save matlab workspace for future reference
save('ica_results.mat')

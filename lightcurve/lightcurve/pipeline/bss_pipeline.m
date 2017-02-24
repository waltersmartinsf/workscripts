%General structure of the programme. This defines the flow of the pipeline
%most functions are contained in sub-modules 



%loading parameters into pipeline parameter stream, PSTREAM
fprintf('loading parameters from config file...')
clearvars -global
global PSTREAM;
PSTREAM = struct();
PSTREAM = bss_configure(PSTREAM);
fprintf('done \n')

%loading data into pipeline data stream, DSTREAM, 
fprintf('Loading data into data-stream...')
global DSTREAM
DSTREAM = struct();
%DSTREAM.DATA.raw = load([PSTREAM.INPUT.path PSTREAM.INPUT.filename]);
DSTREAM.DATA.raw = X;
[Xs1,Xs2] = size(DSTREAM.DATA.raw);
if Xs1 < Xs2 
    DSTREAM.DATA.raw = transpose(DSTREAM.DATA.raw);
end

% DSTREAM.DATA.phase = PAR.PHASE;
DSTREAM.DATA.phase = 0;
fprintf('done \n')





%optional ksr
if strcmp(PSTREAM.KERNEL.onoff,'on') == 1
    fprintf('running kernel smoothing...')
    [DSTREAM,PSTREAM] = bss_ksr(DSTREAM,PSTREAM);
    fprintf('done \n')
end

%running PCA analysis
fprintf('running PCA de- and re-composition algorithm...')
[DSTREAM, PSTREAM] = bss_pca(DSTREAM,PSTREAM);
fprintf('done \n')

%running multi-combi on data
fprintf('running multi-combi algorithm... \n')
[DSTREAM, PSTREAM] = bss_mcombi(DSTREAM,PSTREAM);
fprintf('done \n')

%running signal-identification on separated components
fprintf('running signal identification algorithm...')
[DSTREAM, PSTREAM] = bss_find(DSTREAM,PSTREAM);
fprintf('done \n')

%creating noise model
fprintf('creating noise model...')
[DSTREAM,PSTREAM] = bss_nmodel(DSTREAM,PSTREAM);
fprintf('done \n')


%saving data-stream
% fprintf('saving data-stream to \"DSTREAM.mat\"...')
% save('DSTREAM.mat','DSTREAM')
% fprintf('done \n')



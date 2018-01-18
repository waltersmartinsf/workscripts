% this just wraps the bbs_pipeline code in a function wrapper to be used with pymatbridge

function [components, ISRefica, ISRwasobi, Wmatrix, Wefica, Wwasobi, NoiseModel, nonseparablecomponents, NoiseSignal, NoiseDIAG] = icapy(args)
	%General structure of the programme. This defines the flow of the pipeline
	%most functions are contained in sub-modules

	signals = args %Walter cahnge this line

	%loading parameters into pipeline parameter stream, PSTREAM
	clearvars -global
	global PSTREAM;
	PSTREAM = struct();
	PSTREAM = bss_configure(PSTREAM);

	%loading data into pipeline data stream, DSTREAM,
	%fprintf('Loading data into data-stream...')
	global DSTREAM
	DSTREAM = struct();
	%DSTREAM.DATA.raw = load([PSTREAM.INPUT.path PSTREAM.INPUT.filename]);
	DSTREAM.DATA.raw = signals;
	[Xs1,Xs2] = size(DSTREAM.DATA.raw);

	if Xs1 < Xs2
	    DSTREAM.DATA.raw = transpose(DSTREAM.DATA.raw);
	end

	% DSTREAM.DATA.phase = PAR.PHASE;
	DSTREAM.DATA.phase = 0;

	%optional ksr
	if strcmp(PSTREAM.KERNEL.onoff,'on') == 1
	    %fprintf('running kernel smoothing...')
	    [DSTREAM,PSTREAM] = bss_ksr(DSTREAM,PSTREAM);
	end

	%running PCA analysis
	%fprintf('running PCA de- and re-composition algorithm...')
	[DSTREAM, PSTREAM] = bss_pca(DSTREAM,PSTREAM);

	%running multi-combi on data
	%fprintf('running multi-combi algorithm... \n')
	[DSTREAM, PSTREAM] = bss_mcombi(DSTREAM,PSTREAM);

	%running signal-identification on separated components
	%fprintf('running signal identification algorithm...')
	[DSTREAM, PSTREAM] = bss_find(DSTREAM,PSTREAM);

	%creating noise model
	%fprintf('creating noise model...')
	[DSTREAM,PSTREAM] = bss_nmodel(DSTREAM,PSTREAM);

	components = DSTREAM.DATA.MCOMBI.signals
    nonseparablecomponents = DSTREAM.DATA.MCOMBI.nonsep %edited by walter
    ISRefica = DSTREAM.INFO.MCOMBI.ISRef1 %edited by walter
    ISRwasobi = DSTREAM.INFO.MCOMBI.ISRwa1 %edited by walter
    Wmatrix = DSTREAM.DATA.MCOMBI.w %edited by walter
    Wefica = DSTREAM.INFO.MCOMBI.Wefica %edited by walter
    Wwasobi = DSTREAM.INFO.MCOMBI.Wwasobi %edited by walter
    NoiseModel = DSTREAM.DATA.NMODEL.noise %edited by walter (noise model)
	NoiseSignal = DSTREAM.DATA.NMODEL.signal %edited by walter (noise signal)
	NoiseDIAG = DSTREAM.DATA.NMODEL.nscale %edited by walter (diagonal to scale the noise model)

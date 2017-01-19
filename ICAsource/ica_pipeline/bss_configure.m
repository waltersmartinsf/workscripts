%configuration file for bss pipeline. The PSTREAM structure is defined here

function PSTREAM = bss_configure(PSTREAM)

%general parameters for various modules of the bss_pipeline

PSTREAM.INPUT.path = '/Users/ingowaldmann/Desktop/H_cohe_2';

PSTREAM.INPUT.filename = 'night5.mat';

PSTREAM.PCA.maxcoeff = 0; %sets the number of PCAs retained in data. 0 = no dimensionality reduction, -1 = code decides itself what to do. 

PSTREAM.MCOMBI.arcoeff = 5; %sets max number of AR coefficients for wasobi

PSTREAM.FIND.pearsonthr = 0.7; %detection threshold for lightcurve signal component detection

PSTREAM.KERNEL.onoff = 'off'; %switching kernel smoothing on or off. Default: 'off'

PSTREAM.KERNEL.bandwidth = 0.005; %bandwidth of kernel smoothing 

PSTREAM.LBQ.alpha = 0.05; %lijung-box test alpha value (default: 0.05)

PSTREAM.FIT.lcindex = 1; %lightcurve to be fitted (index of the DSTREAM.DATA.raw matrix)



%parameters for mandel & agol lightcurve fitting only Rp/Rs is fitted 
% 
% PSTREAM.FIT.LCPARAMS.period = 2.21857312; %period of planet (days)
% 
% PSTREAM.FIT.LCPARAMS.epoch = 2453808.91682; %discovery epoch (JD)
% 
% PSTREAM.FIT.LCPARAMS.sma = 0.03142; %orbital semi-major axis (AU)
% 
% PSTREAM.FIT.LCPARAMS.inc = 85.51; %orbital inclination (deg.)
% 
% PSTREAM.FIT.LCPARAMS.ecc = 0; %orbital eccentricity 
% 
% PSTREAM.FIT.LCPARAMS.omega = 0; %argument of periastron (deg)
% 
% PSTREAM.FIT.LCPARAMS.rstar = 0.788; %radius of star (R(solar))
% 
% PSTREAM.FIT.LCPARAMS.u1 = 0.0104; %quadratic limb darkening parameter 1
% 
% PSTREAM.FIT.LCPARAMS.u2 = 0.4663; %quadratic limb darkening parameter 2



PSTREAM.FIT.LCPARAMS.period = 3.9415128; %period of planet (days)

PSTREAM.FIT.LCPARAMS.epoch = 2453808.91682; %discovery epoch (JD)

PSTREAM.FIT.LCPARAMS.sma = 0.0488; %orbital semi-major axis (AU)

PSTREAM.FIT.LCPARAMS.inc = 89.31; %orbital inclination (deg.)

PSTREAM.FIT.LCPARAMS.ecc = 0; %orbital eccentricity 

PSTREAM.FIT.LCPARAMS.omega = 0; %argument of periastron (deg)

PSTREAM.FIT.LCPARAMS.rstar = 0.928; %radius of star (R(solar))

PSTREAM.FIT.LCPARAMS.u1 = 0.05; %quadratic limb darkening parameter 1

PSTREAM.FIT.LCPARAMS.u2 = 0.35; %quadratic limb darkening parameter 2
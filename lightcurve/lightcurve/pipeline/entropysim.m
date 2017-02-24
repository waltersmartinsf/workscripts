%small simulation routine for entropy paper 

function [MODEL2,Gauss, Semigauss] = bss_simu(varargin)

options.default1 = 1;
options.default2 = 1;
for i=1:2:length(varargin)
    switch varargin{i}
        case 'LCPARAM'
            options.LCPARAM = varargin{i+1};
            options.default1 = 0;
        case 'A'
            options.A = varargin{i+1};
            options.default2 = 0;
       otherwise
            warning(['Unknown option: ' varargin{i}])
    end
end

%default planetary parameters for mandel & agol fit 
LCPARAM.Rplanet = 1.151; %GJ1214:0.238917;  %in Rjupiter 
LCPARAM.Rstar = 0.788; %GJ1214: 0.211;     %in Rsun
LCPARAM.INC = 85.67; %GJ1214: 88.62;       %inclination in degrees
LCPARAM.SMA = 0.03142; %GJ1214: 0.01438;      %semi-major-axis in AU
LCPARAM.TIME.STARTPHASE = 0.470;
LCPARAM.TIME.ENDPHASE = 0.530;
LCPARAM.TIME.DATANUM = 1000;

if options.default1 == 0
    LCPARAM = options.LCPARAM;
end


%creating lightcurve phase
LCPARAM.PHASE = zeros(LCPARAM.TIME.DATANUM,1);

step = (LCPARAM.TIME.ENDPHASE - LCPARAM.TIME.STARTPHASE) ./ LCPARAM.TIME.DATANUM;

LCPARAM.PHASE(1,1) = LCPARAM.TIME.STARTPHASE;
for i=2:length(LCPARAM.PHASE)
    LCPARAM.PHASE(i,1) = LCPARAM.PHASE(i-1,1) + step;
end

%creating lightcurve model 
MODEL = agol(LCPARAM.Rplanet,LCPARAM.Rstar,LCPARAM.INC,LCPARAM.SMA,LCPARAM.PHASE(:,1),0,0,0.1);
LCPARAM.MODEL = MODEL;

%%%creating systematic noise now 

%sine curve noise 
sinecur = sin(LCPARAM.PHASE*350.) .* 0.001 + 1.0;

%saw tooth function noise
sawtoothcur = sawtooth(LCPARAM.PHASE*250.) * 0.004 +1.0;

%time-correlated (AR) 
AR=[1 -2.7607 3.8106 -2.6535 0.9238];
% AR(4) coefficients
arcurve=filter(1,AR,0.00004*randn(LCPARAM.TIME.DATANUM,1))+1.0; 

%creating Gaussian noise
Gauss = 0.001 * randn(LCPARAM.TIME.DATANUM,1) + 1.0;
%Gauss2 = 0.001 * randn(length(LCPARAM.PHASE),1) + 1.0;
%Gauss3 = 0.001 * randn(length(LCPARAM.PHASE),1) + 1.0;

MODEL2 = MODEL + sinecur + sawtoothcur + arcurve + Gauss;
MODEL2 = MODEL2 ./ mean(MODEL2);


Semigauss = Gauss + sawtoothcur;
Semigauss = Semigauss ./ mean(Semigauss);


% Module creating simulated data sets for the BSSpipeline 

% Ingo Waldmann  2011

%%%%%%%%%%%%%%%%%%%%%%
%running this module will produce 5 mixed signals using a random or 
%predefined mixing matrix A
%the 5th pre-mixed component will be a three stage AR process testing the
%pipeline for temporarily correlated noise
%
%INPUT:
%LCPARAM = structure containing lightcurve parameters, if not given,
%   HD189733b secondary eclipse parameters are assumed
%A = predefined mixing ratio, if not defined, a random mixing ratio is
%   generated
%
%
%OUTPUT:
%X = structure of mixed signals 
%PREMIX = structure of signals before mixing 
%A = mixing matrix used
%PARAM = updated set of lightcurve and system parameters
%%%%%%%%%%%%%%%%%%%%%%

function [X,PREMIX,A,PARAMOUT] = bss_simu(varargin)

options.default1 = 1;
options.default2 = 1;
options.default3 = 1;
G=0;
for i=1:2:length(varargin)
    switch varargin{i}
        case 'LCPARAM'
            options.LCPARAM = varargin{i+1};
            options.default1 = 0;
        case 'A'
            options.A = varargin{i+1};
            options.default2 = 0;
        case 'G'
            G = varargin{i+1};
            options.default3 = 0;
       otherwise
            warning(['Unknown option: ' varargin{i}])
    end
end

G

%default planetary parameters for mandel & agol fit 
LCPARAM.Rplanet = 1.151; %GJ1214:0.238917;  %in Rjupiter 
LCPARAM.Rstar = 0.788; %GJ1214: 0.211;     %in Rsun
LCPARAM.INC = 85.67; %GJ1214: 88.62;       %inclination in degrees
LCPARAM.SMA = 0.03142; %GJ1214: 0.01438;      %semi-major-axis in AU
LCPARAM.TIME.STARTPHASE = 0.47;%0.470;
LCPARAM.TIME.ENDPHASE = 0.530;
LCPARAM.TIME.DATANUM = 500;

if options.default1 == 0
    LCPARAM = options.LCPARAM;
end


%creating lightcurve phase
LCPARAM.PHASE = zeros(LCPARAM.TIME.DATANUM,1);
ALTPHASE = zeros(LCPARAM.TIME.DATANUM,1);
ALTPHASE(1) = 1;

step = (LCPARAM.TIME.ENDPHASE - LCPARAM.TIME.STARTPHASE) ./ LCPARAM.TIME.DATANUM;

LCPARAM.PHASE(1,1) = LCPARAM.TIME.STARTPHASE;
for i=2:length(LCPARAM.PHASE)
    LCPARAM.PHASE(i,1) = LCPARAM.PHASE(i-1,1) + step;
    ALTPHASE(i) = i;
end

%creating lightcurve model 
LCPARAM.MODEL = agol(LCPARAM.Rplanet,LCPARAM.Rstar,LCPARAM.INC,LCPARAM.SMA,LCPARAM.PHASE(:,1),0,0,0.02);


%%%creating systematic noise now 

%sine curve noise 
sinecur = sin(LCPARAM.PHASE*600.) .* 0.002 + 1.0;

%saw tooth function noise
sawtoothcur = sawtooth(LCPARAM.PHASE*250.) * 0.004 +1.0;

%time-correlated (AR) 
AR=[1 -2.7607 3.8106 -2.6535 0.9238];
% AR(4) coefficients
arcurve=filter(1,AR,0.00004*randn(LCPARAM.TIME.DATANUM,1))+1.0; 

%log curve
c = [0.005,3,100,4,5];
expcurv = c(1).*(1 - c(2).*exp(-ALTPHASE/c(3)) - c(4).*exp(-ALTPHASE/c(5)));


%creating Gaussian noise
Gauss = G * randn(length(LCPARAM.PHASE),1) + 1.0;
% Gauss2 = G * randn(length(LCPARAM.PHASE),1) + 1.0;
%Gauss3 = 0.001 * randn(length(LCPARAM.PHASE),1) + 1.0;


%adding gaussian noise to individual channels
% if options.default3 == 0
    Gauss2 = 0.0001 .* randn(length(LCPARAM.PHASE),1);
    LCPARAM.MODEL = LCPARAM.MODEL + Gauss2;
    sinecur = sinecur + Gauss2;
    sawtoothcur = sawtoothcur + Gauss2;
% end


%creating umixed data matrix
S = cat(2,LCPARAM.MODEL,sinecur,sawtoothcur,arcurve,Gauss);
%S = cat(2,sinecur,sawtoothcur,arcurve,Gauss);

%loading signals into PREMIX structure.
PREMIX.lc = LCPARAM.MODEL;
PREMIX.sine = sinecur;
PREMIX.saw =sawtoothcur;
PREMIX.ar = arcurve;
% PREMIX.expcurv = expcurv;
PREMIX.gauss = Gauss;
PREMIX.matrix = S;

%creating random mixing matrix
if options.default2 == 0
    A = options.A;
else
    %A = rand(4,4);
    A = rand(5,5);
end
A = A * 10;
save('Amix.mat','A')

%mixing the signals X = AS
X = transpose(A*S');

[s1,s2] = size(X);

for i=1:s2
    %X(:,i) = X(:,i) + LCPARAM.MODEL;
    X(:,i) = X(:,i) ./ mean(X(:,i));
end

%filling PARAMOUT structure
PARAMOUT = LCPARAM;

%plotting generated signals

% figure(3)
% plot(Gauss)
% hold on
% plot(LCPARAM.MODEL)

% figure(1)
% subplot(5,1,1)
% plot(LCPARAM.PHASE,LCPARAM.MODEL)
% title('individual unmixed components')
% subplot(5,1,2)
% plot(LCPARAM.PHASE,sinecur)
% subplot(5,1,3)
% plot(LCPARAM.PHASE,sawtoothcur)
% subplot(5,1,4)
% plot(LCPARAM.PHASE,arcurve)
% subplot(5,1,5)
% plot(LCPARAM.PHASE,Gauss)
% xlabel('phase')
% ylabel('simulated relative flux')
% hold off
% 
% 
% figure(2)
% plot(LCPARAM.PHASE, X(:,1),'x')
% hold on
% plot(LCPARAM.PHASE, X(:,2),'rx')
% plot(LCPARAM.PHASE, X(:,3),'gx')
% plot(LCPARAM.PHASE, X(:,4),'kx')
% plot(LCPARAM.PHASE, X(:,5),'yx')
% title('mixed components')
% xlabel('phase')
% ylabel('simulated relative flux')
% hold off
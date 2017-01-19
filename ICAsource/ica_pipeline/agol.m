% Matlab implementation of mandel and agol curve for secondary eclipses
%ie small planet approximation 

% ;+
% ; PURPOSE: make a planet light curve according to the article of
% ; Mandel & Agol 2002, ApJ 580
% 
% ; INPUT:
% ;  - radplanet:    radius of the planet [cm]
% ;  - radstar:      radius of the star   [cm]
% ;  - inc:          the inclination [rad]
% ;  - sma:          semimajor axis [cm]
% ;  - phase:        the phase to calculate the curve at [NORMALIZED]
% ; OPTIONAL input:
% ;  - eclipsedepth: what is the real depth of the eclipse to be modelled
% ;- 

%%%%%%%%%%%%%%%%%

function [FluxRatio] = agol(radplanet,radstar,inclination,semimajoraxis,phase,ecc,omega,depth)

%converting units to cgs

inclination = inclination / 360 * 2*pi;
semimajoraxis = semimajoraxis *149.6e10;
radstar = radstar * 6.955e9;
radplanet = radplanet * 7.1492e8;



% go from phase to center-to-center distance
%distancevector = delta((phase-0.5d0)*2d0*!dpi,inclination)*semimajoraxis/radstar

distancevector = delta(((phase)*2d0*pi),inclination,ecc,omega)*semimajoraxis/radstar;


% variables defined in Paper Mandel & Agol 2002

  rp   = double(radplanet);
  rs   = double(radstar);
  z    = double(distancevector);
  p    = double(rp/rs);
  
  %inclination = double(inclination);
  %semimajoraxis = double(semimajoraxis);
  
%;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
%;;
%;; The planet light curve for a uniform source
%;;
%;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


lambda = zeros(size(z));
kappa_1 = acos((1-p.^2+z.^2)./2./z);
kappa_0 = acos((p.^2+z.^2-1)./2./p./z);


idx = find(z > (1 + p));
if (~isempty(idx)),
    lambda(idx) = 0;
end

idx = find(z > abs(1-p));
if (~isempty(idx)),
    lambda(idx,1) = 1 ./ pi * (p.^2.*kappa_0(idx,1)+kappa_1(idx,1)- ...
        sqrt((4*z(idx,1).^2-(1+z(idx,1).^2 - p.^2).^2)/4));
end


idx = find(z < (1+p));
if (~isempty(idx)),
    lambda(idx) = 1 ./ pi * (p.^2.*kappa_0(idx)+kappa_1(idx)- ...
        sqrt((4*z(idx).^2-(1+z(idx).^2-p.^2).^2)/4));
end


idx = find(z < (1-p));
if (~isempty(idx)),
    lambda(idx) = p.^2;
end

idx = find(z < (p-1));
if (~isempty(idx)),
    lambda(idx) = 1;
end

FluxRatio = 1 - lambda;

if (depth ~= 0),
    scale = depth ./ (rp./rs).^2;
    lambda = lambda * scale;
    FluxRatio = 1 - lambda;
end

FluxRatio = real(FluxRatio);






















% ;+
% ; PURPOSE: compute the center-to-center distance between star and
% ; planet given the orbital parameters. 
% ;
% ; INPUT:
% ;  - theta = orbital phase [ rad ! ]
% ;  - inc   = inclination of the system [rad]
% ; OPTIONAL INPUT:
% ;  - ecc   = eccentricity (default = 0)
% ;  - omega = omega if the eccentricity ne 0.
% ;- 

function [delta] = delta(theta,inc,ecc,omega)

%poo = (1 - cos(theta)) ..^ 2

delta = sqrt(1-((cos(theta)).^2*(sin(inc)).^2));

% delta = sqrt((sin(theta)).^2 + (cos(theta).*cos(theta)).^2);


if (ecc ~= 0) || (omega ~= 0),
    thetan = theta;
    delta = (1-ecc.^2)/(1-ecc*sin(thetan-omega))* ...
        sqrt((1-(cos(thetan)).^2*(sin(inc)).^2));
end


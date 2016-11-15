function [r,theta,phi] = cart2sph(x,y,z)
% Function to convert coordinates from Cartesian to Spherical system.
%
% TJ Young
% 14 November 2016

%% Conversion to spherical coordinates
r = sqrt(x.^2 + y.^2 + z.^2); % Radius
theta = acosd(abs(z)./r); % Polar angle
phi = atand(y./x); % Azimuthal angle

%% Correct ambiguous values of phi (as period of tangent is pi)
for ii = 1:numel(phi)
    % Correct values located in 2nd quadrant
    if x(ii) < 0 && y(ii) > 0
        phi(ii) = phi(ii) + 180; 
    % Correct values located in 3rd quadrant
    elseif x(ii) < 0 && y(ii) < 0 
        phi(ii) = phi(ii) + 180; 
    % Correct values located in 4th quadrant
    elseif phi(ii) < 0
        phi(ii) = phi(ii) + 360; 
    end
end
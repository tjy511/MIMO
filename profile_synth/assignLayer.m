function [zz,zz_reflect] = assignLayer(xx,yy,zz,zz_reflect,L,A,bedSelect,param,addBed)

% Implements a layer in a 2D ice profile through depth and assigns a
% reflectivity parameter to all identified pixels
%
% TJ Young
% 09 December 2016



switch bedSelect
    case 'flat'
        param = []; % Dummy parameter
        idx = find(yy(:,1) == L); % Find pixels at reflector depth
        zz(idx,:) = db2mag(A); % Change value of pixels to layer reflectivity
        zz_reflect(idx,:) = 0; % Assigns flat reflector (0 degrees) to ref array
        
    case 'slope'
        md = param; % Slope in degrees
        m = tand(md); % Convert to slope format 
        lSlope = m .* xx(1,:) + L; % Assigns sloped reflector
        [~,idx] = min(abs(bsxfun(@minus,yy,lSlope))); % Find nearest pixels to reflector depth
        for ii = 1:length(idx)
            zz(idx(ii),ii) = db2mag(A); %  Change value of pixels to layer reflectivity
            zz_reflect(idx(ii),ii) = -md; % Assigns slope value to ref array
        end
        
    case 'sine'
        amp = param(1); period = param(2); % Parameters
        lSin = amp.*sind(xx(1,:) .* period .* 360/size(xx,2)) + L; % Assigns sine reflector
        [~,idx] = min(abs(bsxfun(@minus,yy,lSin))); % Find nearest pixels to reflector depth
        grad = gradient(lSin); gradtan = -atand(grad); % Finds slope of layer (degrees); note that depth is positive and thus a negative is neccesary
        for ii = 1:length(idx)
            zz(idx(ii),ii) = db2mag(A); % Change value of pixels to layer reflectivity
            zz_reflect(idx(ii),ii) = gradtan(ii); % Assigns reflector array to ref array
        end
end

% Add bed
if exist('addBed','var')
    Asg = A-20; % Power of subglacial profile
    for ii = 1:length(idx)
        zz(idx(ii)+1:end,ii) = db2mag(Asg); % Change value of pixels to subglacial layer reflectivity
        zz_reflect(idx(ii)+1:end,ii) = 0; % Assigns flat reflector (0 degrees) to ref array
        
        noise = randn(size(zz(idx(ii)+1:end,ii)));
        noise = 100 .* noise .* db2mag(Asg);
        zz(idx(ii)+1:end,ii) = zz(idx(ii)+1:end,ii)+noise;
        
    end
end
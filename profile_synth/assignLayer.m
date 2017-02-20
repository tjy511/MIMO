function [zz,zz_reflect] = assignLayer(xx,yy,zz,zz_reflect,L,A,select,param,varargin)

% Implements a layer in a 2D ice profile through depth and assigns a
% reflectivity parameter to all identified pixels
%
% TJ Young
% 09 December 2016

if nargin < 9
    bedAdd1 = 'none'; bedAdd2 = 'none';
end

if nargin >= 9
    bedAdd1 = varargin{1};
    bedAdd2 = 'none';
    if nargin == 11
        bedAdd2 = varargin{2};
        bedParam = varargin{3};
    end
end

switch select
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

% Add subglacial mountain (unfinished)
if strcmp(bedAdd1,'addHill') | strcmp(bedAdd2,'addHill')
    try
        location = bedParam(1); % Location of mountain centre
        wid = bedParam(2); % Width of mountain
        apex = bedParam(3); % Height above bed
        % Find indices of mountain
        idxx = ceil(size(xx,2)*(location-wid/2)) : ceil(size(xx,2)*(location+wid/2)); % x index for mountain
        idxAx = floor(length(idxx)/2) + idxx(1); % x index for mountain apex
        [~,idxAy] = min(abs(yy(:,1)-(L-apex))); % y index for mountain apex
        idxy = interp1([idxx(1) idxAx idxx(end)],[idx(idxx(1)) idxAy idx(idxx(end))],idxx); % Interpolate between points to obtain mountain
        idxy = round(idxy); % Round to obtain y index for mountain
        % Find gradient of mountain
        gradm = (idxy(1)-idxAy) / (idxx(1)-idxAx); % Slope of mountain (1 side)
        gradm = [repmat(gradm,[1 floor(length(idxx)/2)]) 0 -repmat(gradm,[1 floor(length(idxx)/2)])]; % Extrude slope to indices
        gradmtan = -atan(gradm); % Finds slope of layer (degrees); note that depth is positive and thus a negative is neccesary
        % Rewrite pixels with identified indices
        for ii = 1:length(idxx)
            zz(idxy(ii),idxx(ii)) = db2mag(A); % Change value of pixels to layer reflectivity
            zz_reflect(idxy(ii),idxx(ii)) = gradmtan(ii); % Assigns reflector array to ref array
        end
        idx(idxx) = idxy; % Rewrite original index to include mountain
    end
end

% Add subglacial bed
if strcmp(bedAdd1,'addBed') | strcmp(bedAdd2,'addBed')
    try
        Asg = A-20; % Power of subglacial profile
        for ii = 1:length(idx)
            zz(idx(ii)+1:end,ii) = db2mag(Asg); % Change value of pixels to subglacial layer reflectivity
            for jj = idx(ii)+1:size(zz_reflect,1)
                zz_reflect(jj,ii) = rand(1)*180-90; % Assigns random reflector to ref array
            end
            %zz_reflect(idx(ii)+1:end,ii) = 0; % Assigns zero reflector to ref array
            
            noise = randn(size(zz(idx(ii)+1:end,ii)));
            noise = 100 .* noise .* db2mag(Asg);
            zz(idx(ii)+1:end,ii) = zz(idx(ii)+1:end,ii)+noise;
        end
    end
end
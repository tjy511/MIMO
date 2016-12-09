function [zz,zz_reflect] = assignBed(xx,yy,zz,zz_reflect,L,A,bedSelect)

% Implements a bed in a 2D ice profile through depth and assigns a
% reflectivity parameter to all identified pixels
%
% TJ Young
% 09 December 2016

%zzOrig = zz; % TEMP
%bedSelect = 'sine';
%%

%zz = zzOrig;
switch bedSelect
    case 'flat'
        idx = find(yy(:,1) == L); % Find pixels at reflector depth
        zz(idx,:) = db2mag(A); % Changes value of pixels to bed reflectivity
        zz_reflect(idx,:) = 0; % Assigns flat reflector (0 degrees) to ref array
    case 'sine'
        amp = pi; period = 10; 
        bedSin = pi.*sind(xx(1,:) .* period .* 360/size(xx,2)) + L; 
        idx = find(xx(1,:) == bedSin); % Need ot use nearest function
end


%% SCRATCH
%figProfileSynth = plotimgprofile_gland(xx,yy,zz);
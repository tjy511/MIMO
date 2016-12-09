% Script to test sensitivities of parameters to 2D power return
%
% Current: Isotropic bed reflector
%
%
% TJ Young
% 07 December 2016

ccc

%% Parameters

% Depth profile
N = -100; % Noise floor level (dB)
A = -20; % Maximum reflector reflectivity (dB)
L = 600; % Depth of reflector (m)

% Reflectivity
sigma = (68.27/2)  / (10/9); % Standard deviation (bell curve) [degrees]

% Bed formation
bedSelect = 'flat';

% Beam forming
thetaStA = -30:30; % Steering angle range
antSelect = 'pencil'; % Antenna selection
freq = 3e8; % MIMO frequency (in 'bowtie')

%% Synthesise ice profile through depth

% Dimensions of profile
depths = 0:700; widths = -300:300;
xx = repmat(widths,length(depths),1);
yy = repmat(depths',1,length(widths));
zz = repmat(db2mag(N),size(xx));

% Add noise
if ~isnan(N)
    noise = randn(size(zz));
    noise = 100 .* noise .* db2mag(N);
    zz = zz+noise;
end

%% Create pixel distribution

% Create Gaussian distribution of power return
px = -90:1:90; % Slope angle
py = normpdf(px,0,sigma);
py_max = max(py);
py = py .* (1/py_max); % Power return scaled to dB

% Plotting fancies
figPixPowerDist = figure;
plot(px,py)
xlim([px(1) px(end)])
xlabel(['Slope angle [' 176 ']'])
ylabel('Power return [dB]')

%% Add layer to ice profile

zz_reflect = nan(size(zz)); % Reflectivity of pixels
[zz,zz_reflect] = assignBed(xx,yy,zz,zz_reflect,L,A,'flat'); 
figProfileSynth = plotimgprofile_gland(xx,yy,zz);

%% Sum through all ranges

R_loc = [1 ceil(length(widths)/2)]; % Location of radar on zz
zz_r = sqrt(xx.^2 + yy.^2); % Distance of pixels from radar
zz_theta = asind(xx./zz_r); % Angle of pixels from radar

% Find cells through all depths
clear R_sum R_sum_norm
R_sum = nan(length(thetaStA),length(depths));

% Loop through steering angles
for ii = 1:length(thetaStA)
    % Load antenna radiation pattern
    [thetaScA,~,antRE] = antennaBP(antSelect,0,thetaStA(ii));
    
    % Load array factor
    switch antSelect
        case 'pencil'
            arrRE = antRE;
        case 'bowtie'
            pos = array2d_final_gland(0); % Calculate virtual antenna locations
            w = ones(1,size(pos,1))/size(pos,1); % Calculate weights
            c = 3e8/sqrt(3.18); 
            W = arrayFactor(pos(:,1),pos(:,2),w,freq,c,thetaScA,0,thetaStA(ii)); % Unweighted
            arrRE = W' .* antRE; % Weight array factor with antenna radiation pattern
    end
    
    % Loop through depths
    for jj = 1:length(depths)
        clear R_idx R_pix R_w
        
        % Locate pixels at specified range from source
        depth = zz_r(jj,R_loc(2)); % Depth vector
        R_idx = findRangePixels(zz_r,depth); % Finds pixels at specified depth
        R_pix(:,1) = zz(R_idx); % Power of identified pixels
        R_pix(:,2) = zz_r(R_idx); % Distance of identified pixels from point
        R_pix(:,3) = zz_theta(R_idx); % Angle of identified pixels from point
        R_pix(:,4) = zz_reflect(R_idx); % Identification of layers in identified pixels
        
        clear w_rad_idx R_w
        wrad_idx = nan(size(R_pix,1),1); R_wrad = wrad_idx;
        R_wref = ones(size(R_pix,1),1);
        for kk = 1:size(R_pix,1)
            
            % Weight identified pixels to antenna radiation pattern
            [~,wrad_idx(kk,1)] = min(abs(thetaScA-R_pix(kk,3)));
            R_wrad(kk,1) = arrRE(wrad_idx(kk,1));
            
            % Weight identified pixels to layer reflectivity
            if ~isnan(R_pix(kk,4))
                [~,wref_idx] = min(abs(px-R_pix(kk,3)));
                R_wref(kk,1) = py(wref_idx);
            end
        end
        
        % Sum all pixels at specified range from source
        R_sum(ii,jj) = sum(R_pix(:,1) .* R_wrad(:,1).^2 .* R_wref(:,1));
    end
    
    disp(['Processed power return at steering angle: ',num2str(thetaStA(ii)),' degrees'])
end

%% Plotting fancies

% Plot profile return distribution (1D)
figure
plot(depths,db(R_sum(ceil(length(thetaStA)/2),:)))
ylim([-100 20])
xlabel('Depth [m]')
ylabel('Power return [dB]')

% Set up 2D grid for profile distribution
[rr,tt] = meshgrid(-depths,linspace(min(thetaStA),max(thetaStA),length(thetaStA)));
rx = rr.*sind(tt);
ry = rr.*cosd(tt);
rz = R_sum;

% Plot profile return distribution (2D)
figProfileReturn = plotimgdepth_gland(rx,ry,db(rz));
axis equal
caxis([-120 20])
legend = colorbar('Ticks',[-120 -100 -80 -60 -40 -20 0 20]);

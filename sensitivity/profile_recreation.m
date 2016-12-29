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
N = -100; % Noise floor level [dB]
Al = -40; % Internal reflector reflectivity [dB]
Ab = -20; % Bed reflector reflectivity [dB]
L = 300; % Depth of internal reflector [m]
B = 600; % Depth of bed reflector [m]

% Reflectivity
sigma = (68.27/2)  / (10/9); % Standard deviation (bell curve) [degrees]

% Bed formation
bedSelect = 'sine'; % 'flat' 'slope' 'sine'

% Antenna location
antLoc = 'other'; % 'store1' 'other'
txrx = [32 1]; % Default: [8 8]
quadrant = 4; % Default: 4
off = [2.49 -2.49]; % Default: [2.49 -2.49] [m]
dPhy = 0.83; % Default: 0.83 [m]

% Beam forming
freq = 3e8; % MIMO frequency [Hz]
thetaStA = -30:30; % Steering angle range [degrees]
antSelect = 'pencil'; % Antenna selection ('isotropic' 'pencil' 'bowtie' 'helix' 'dipole')
weight = 1; % Weighting beam with array factor [binary]

% Switches
doPlot = 1; 
doSave = 1;

%% Create pixel distribution

% Create Gaussian distribution of power return
px = -90:1:90; % Slope angle
py = normpdf(px,0,sigma);
py_max = max(py);
py = py .* (1/py_max); % Power return scaled to dB

% Plotting fancies
fig_pwr = figure;
plot(px,py)
xlim([px(1) px(end)])
xlabel(['Slope angle [' 176 ']'])
ylabel('Power return [dB]')

%% Synthesise ice profile through depth

% Dimensions of profile
depths = 0:1:700; widths = -300:1:300;
xx = repmat(widths,length(depths),1);
yy = repmat(depths',1,length(widths));
zz = repmat(db2mag(N),size(xx));
zz_reflect = nan(size(zz)); % Reflectivity of pixels

% Add internal layers to profile
[zz,zz_reflect] = assignLayer(xx,yy,zz,zz_reflect,100,Al,'flat');
[zz,zz_reflect] = assignLayer(xx,yy,zz,zz_reflect,200,Al,'slope',2);
[zz,zz_reflect] = assignLayer(xx,yy,zz,zz_reflect,300,Al,'slope',5);
[zz,zz_reflect] = assignLayer(xx,yy,zz,zz_reflect,400,Al,'slope',10);

% Add bed layer to profile
[zz,zz_reflect] = assignLayer(xx,yy,zz,zz_reflect,B,Ab,'sine',[pi 10],1); 

% Add noise to profile
if ~isnan(N)
    noise = randn(size(zz));
    noise = 100 .* noise .* db2mag(N);
    zz = zz+noise;
end

% Plot profile
fig_profsynth = plotimgprofile_gland(xx,yy,zz);

%% Sum through all ranges

% Parameters of zz
R_loc = [1 ceil(length(widths)/2)]; % Location of radar on zz
zz_r = sqrt(xx.^2 + yy.^2); % Distance of pixels from radar
zz_theta = asind(xx./zz_r); % Angle of pixels from radar

% Calculate virtual antenna locations
if strcmp(antLoc,'store1') == 1
    ve = antennaLoc(antSelect,doPlot);
else
    ve = antennaLoc(txrx,quadrant,off,dPhy,doPlot);
end
xPos = ve(:,1); yPos = ve(:,2);

w = ones(1,size(ve,1))/size(ve,1); % Calculate weights
c = 3e8/sqrt(3.18);

% Find cells through all depths
clear R_sum R_sum_norm
R_sum = nan(length(thetaStA),length(depths));

% Loop through steering angles
for ii = 1:length(thetaStA)
    
    % Load antenna radiation pattern
    [thetaScA,~,antRE] = antennaBP(antSelect,0,thetaStA(ii));
    
    % Load array factor
    W = arrayFactor(xPos,yPos,w,freq,c,thetaScA,0,thetaStA(ii)); % Unweighted
    switch weight
        case 0
            arrRE = antRE; % Unweighted array factor
        case 1
            arrRE = W' .* antRE; % Weight array factor with antenna radiation pattern
    end
    
    % Loop through depths
    for jj = 1:length(depths)
        clear R_pix
        %R_idx = nan(size(zz));
        
        % Locate pixels at specified range from source
        depth = zz_r(jj,R_loc(2)); % Depth vector
        R_idx = findRangePixels(zz_r,depth); % Finds pixels at specified depth
        R_pix = nan(sum(sum(R_idx)),4); % Pre-allocate array
        R_pix(:,1) = zz(R_idx); % Power of identified pixels
        R_pix(:,2) = zz_r(R_idx); % Distance of identified pixels from point
        R_pix(:,3) = zz_theta(R_idx); % Angle of identified pixels from point
        R_pix(:,4) = zz_reflect(R_idx); % Identification of layers in identified pixels
        
        %clear w_rad_idx R_w
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
fig_prof1d = figure;
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
fig_prof2d = plotimgdepth_gland(rx,ry,db(rz));
axis equal
caxis([-120 20])
legend = colorbar('Ticks',[-120 -100 -80 -60 -40 -20 0 20]);

%% Saving fancies
% prof_synth.png prof_return_2d.png 
%set(gcf,'color','w');
%export_fig prof_return_2d.png -m2
%% Saving fancies

if doSave
    % Create and cd to folder
    fileLoc = strcat('~/Google Drive/Academic/papers/paper3/figs/recreation/',...
        num2str(txrx(1)),'-',num2str(txrx(2)),'_',num2str(dPhy*1000));
    try
        cd(fileLoc);
    catch
        mkdir(fileLoc); cd(fileLoc);
    end
    set([fig_prof1d fig_prof2d],'color','w')
    export_fig(fig_prof1d,strcat(antSelect,'_PR_1d.png'),'-m2');
    export_fig(fig_prof2d,strcat(antSelect,'_PR_2d.png'),'-m2');
    if doPlot
        set([fig_pwr fig_profsynth],'color','w')
        export_fig(fig_pwr,'PD_pixel.png','-m2')
        export_fig(fig_profsynth,'PS.png','-m2');
    end
end
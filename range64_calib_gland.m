% Script to account for cable length on MIMO array.
%
% Script written by Lai Bun Lok

%% Set up parameters

% Radar parameters
fc=300e6; % Centre frequency (Hz)
B=200e6; % Bandwidth (Hz)
T=1; % Burst length (s)
c=3e8; % Speed of light (m/s)
er=3.1; % Permittivity of ice (relative to free space)
deltaT=0;
fs=40000; % Samples per subburst
phi0=0; % Offset angle
padfactor=2; % Pad factor (for zero-padding)
k=2*pi*B/T; % Chirp gradient (rad/s^2)
N=64; % total number of virtual elements to consider
Nve=1; % virtual element index to consider
freq=0:1/padfactor:fs/2-(1/padfactor);
rangemetres = freq*c/2/B/sqrt(er); % Targe range
maxrange=1000; % maximum scale on range profile plot

deltaR=c/(2*B*sqrt(er)); % Coarse range bin (m)

% Pre-allocate matrices
chirp=zeros(N,fs); % Deramp data matrix, T1-R1, T1-R2, T1-R3, etc... T8-R8
res=zeros(1,fs);
vtd=zeros(N,fs);
ampspec=zeros(N,fs); % 64xfs matrix to hold all 64 amplitudes of range profiles

% Parameters for antenna cable length calibration
cablelength1=5; % physical length (m) of cable 1
cablelength2=10; % physical length (m) of cable 2
er_cable=2.2; % dielectric constant of the antenna cables

iceequiv_cablelength1=cablelength1*sqrt(er_cable/er);
iceequiv_cablelength2=cablelength2*sqrt(er_cable/er);

rangebinshift1=round(iceequiv_cablelength1*padfactor/deltaR/2);
rangebinshift2=round(iceequiv_cablelength2*padfactor/deltaR/2);

phaseadj1=2*pi*iceequiv_cablelength1*fc/c;
phaseadj2=2*pi*iceequiv_cablelength2*fc/c;
phaseadj1=-1*(rem(phaseadj1+pi,2*pi)-pi); %express in range -pi to pi
phaseadj2=-1*(rem(phaseadj2+pi,2*pi)-pi); %express in range -pi to pi
phaseadj1deg=phaseadj1*180/pi; % Conversion from radians to degrees
phaseadj2deg=phaseadj2*180/pi;

%% DEBUGGING AREA
t=1;
bprev=1;
ind_track=[2420 2420 2421 2413 2416 2425 2426 2426 ...
    2421 2421 2420 2419 2419 2425 2428 2427 ...
    2417 2417 2416 2414 2414 2423 2423 2423 ...
    2418 2417 2417 2430 2414 2423 2423 2423 ...
    2415 2415 2422 2446 2414 2423 2423 2423 ...
    2454 2448 2446 2438 2447 2435 2437 2437 ...
    2445 2451 2438 2439 2445 2443 2443 2442 ...
    2436 2431 2448 2436 2441 2449 2445 2443 ]; % index of the echo to track, 'leading edges' at t=0
b=ind_track;

%% Enter cable configuration for radar
CL_Tx=[5 5 5 5 10 10 10 10]; % Tx1:8 antenna cable lengths (m)
CL_Rx=[5 5 10 10 10 10 10 10]; % Rx1:8 antenna cable lengths (m)

% Work out index of various cable combinations
idx=zeros(64,1);
g=1;
for s=1:8 %Tx
    for q=1:8 %Rx
        if CL_Rx(q)==CL_Tx(s) && CL_Rx(q)<max(CL_Rx)
            idx(g)=1;
        elseif CL_Rx(q)>CL_Tx(s) || CL_Rx(q)<CL_Tx(s)
            idx(g)=2;
        else
            idx(g)=3;
        end
        g=g+1;
    end
end

%% Run through each burst
for kk=1:N
    
    % Sort the data from each virtual element into the chirp array [64x40000] at time index t
    for k=1:64
        chirp(k,1:fs)=vdat.v(1+fs*(k-1):fs*k);
    end
    vif=chirp(kk,:); % select one of the 64 chirp Tx-Rx pair combination to process
    
    % Then apply Brennan's algorithm for each virtual element, at time index t
    for n=1:1
        vtd(n,:) = vif; % store the time domain signal to check for signs of clipping
        vif = vif - mean(vif);
        L = length(vif);
        vif = double(vif).*transpose(blackman(L));
        
        % Zero-pad and shift centre of deramped signal to t=0
        vifpad = zeros(1,padfactor*L);
        vifpad(1:0.5*L) = vif(0.5*L + 1:L);
        vifpad(padfactor*L - 0.5*L + 1:padfactor*L) = vif(1:0.5*L);
        vif = vifpad;
        
        % Take FFT and scale
        fftvif=(sqrt(2*padfactor)/length(vif))*fft(vif);
        spec = fftvif(1:round(length(vif)/2-0.5));
        
        % Apply Brennan's antenna cable calibration procedure to the relevant virtual element combinations   
        if idx(kk)==3 % cable2+cable2
            spec=spec*exp(j*phaseadj2); 
            spec(1:L-rangebinshift2)=spec(1+rangebinshift2:L);
        elseif idx(kk)==2 % cable1+cable2
            spec=spec*exp(j*phaseadj1); 
            spec(1:L-rangebinshift1)=spec(1+rangebinshift1:L);
        else %cable1+cable1
        end
        magspec = abs(spec);
        ampspec(kk,:) = spec;
        
        % Compensate for window
        magspec = magspec/sqrt(mean(blackman(L).*blackman(L)));
        res(kk,:) = magspec; % matrix to hold max 64 FFTs, range profiles
        
        b=ind_track(kk);
        rangebinmetres=(b-1)*c/2/abs(B)/padfactor/sqrt(er);
        
        % Generate reference phase - i.e. expected phase of signal at each range bin
        m=0:length(spec)-1;
        m=m/padfactor;
        phiref=2*pi*fc.*m/B-k*m.*m/2/B/B+phi0;
        deltaphi=2*pi*deltaT*fs*m/length(spec); % Compensate phase for ADC delay of deltaT
        phiref=phiref-deltaphi;
        comp=exp(-j*phiref);
        %spec=spec.*comp; % Compensate measured signal so that phase = 0 in each range bin for target centred in range bin
        
        % Use measured phase for fine-precision range measurement
        phasepeakdegs=(180/pi)*angle(spec(b));
        phasecheck(kk,1)=(180/pi)*angle(spec(2552));
        phasecheck(kk,1)=(180/pi)*angle(spec(2668));
        phasepeakdegscurrent=phasepeakdegs;
        
        % Tracking the phase to +-135deg between range bins
        b=b-1;
        rangebinmetres=(b-1)*c/2/abs(B)/padfactor/sqrt(er);
        phasepeakdegs=(180/pi)*angle(spec(b));
        bprev=b;
        
        deltarangemetres=phasepeakdegs*c/360/2/fc/sqrt(er);
        estimatedrangemetres=rangebinmetres+deltarangemetres;
        
    end
    range_est(1,kk)=estimatedrangemetres;
    
    % Plot the range profile curve from one virtual element point
    % hold all
    % plot(rangemetres,20*log10(res(kk,:))+30-10*log(50),rangebinmetres(t),-80,'+')
    % xlabel('Range (m)')
    % ylabel('Power (dB)')
    % grid
    % %xlim([500 580]) %site 1
    % xlim([400 600]) %site 3
    % ylim([-120 0])
    
    % Plot the estimated basal range from phase-sensitive processing
    % plot(estimatedrangemetres,'+');
    % ylim([498 500])
    
    % Plot the time domain signal from one virtual element node
    %  plot(vtd(1,:));
    %  ylim([0 2.5])
    %  xlim([19900 20100])
    %res_store(t,:)=res;
    %title([datestr(dateStamp(t,1)), ' Virtual element', num2str(kk), ' t=', num2str(t)]);
    %title(['Virtual element', num2str(kk), ' t=', num2str(t)]);
    % M(t)=getframe(gcf);
    t=t+1;
end

% Output data
ampspec_sum=zeros(1,fs);
for kk=1:64
    ampspec_sum=ampspec_sum+ampspec(kk,:);
end
magspecsum=abs(ampspec_sum);
magspecsum = magspecsum/sqrt(mean(blackman(L).*blackman(L)));

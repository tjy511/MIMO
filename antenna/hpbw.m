function hp_bw = hpbw(thetaScA,W)
%
% Calculates half power beam width of radiation pattern. Assumes that the
% main beam is centred at 0 degrees.
%
% Output: Half-power beamwidth
% Inputs: Array factor in dB (W);
%         Theta scanning angle in degrees (thetaScA);
%         Inputs should ideally be centred about 0 degrees.
%
% TJ Young
% 06 December 2016
% Edited 06 January 2016

% Reshape W
if size(W,1) > size(W,2)
    W = W';
end

% Convert theta to degrees if necessary
if range(thetaScA) <= 2*pi
    thetaScA = rad2deg(thetaScA);
end

% Interpolate to improve accuracy
dt = median(diff(thetaScA));
interpTheta = interp1(thetaScA,thetaScA,min(thetaScA):dt/100:max(thetaScA));
for ii = 1:size(W,1)
    WInterp(ii,:) = interp1(thetaScA,W(ii,:),interpTheta);
end

W = WInterp; % Replace W with interpolated W
thetaScA = interpTheta; % Replace thetaScA with interpolated thetaScA
dt = median(diff(thetaScA)); % Update dt

% % Remove duplicate value (0 = 360 degrees)
% if range(thetaScA) == 360
%     W = W(1:end-1);
%     thetaScA = thetaScA(1:end-1);
% end

% Shift main beam to centre of array
warning('off','MATLAB:circshift:ScalarShift')
[~,centreIdx] = min(abs(thetaScA)); % Index of 0 degrees
centreArr = ceil(length(W)/2); % Array centre
shift = centreIdx-centreArr; % # elements to circshift (off-centre if length(W) is even)
W = circshift(W,[0 -shift]); % Shift peak to centre moving right
thetaScA = circshift(thetaScA,[0 -shift]); % Shift 0 degrees to centre moving right

% % Reshape theta
% thetaScA(thetaScA>180 & thetaScA<360) = thetaScA(thetaScA>180 & thetaScA<360)-360; % Reshape theta to be increasingly valued

thresh = 0.5; % HPBW [dB]
WInv = fliplr(W);

% Find lower bound
for ii = 1:size(W,1)
    [~,lb] = find(WInv(ii,centreArr:end) < thresh,1,'first'); % Lower bound index relative to centre
    [~,ub] = find(W(ii,centreArr:end) < thresh,1,'first'); % Upper bound relative to centre
    lb = centreArr-lb; ub = centreArr+ub; % Absolute index
    try
        hp_bw(ii) = (ub-lb) .* dt ./ 2; % HPBW [+/- degrees]
    catch
        hp_bw(ii) = NaN;
    end
end
    
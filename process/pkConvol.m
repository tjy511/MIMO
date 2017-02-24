function vvOut = pkConvol(vvIn,ftype,param)
% 2D Peak Convolution to smooth output radar imagery
% param = [filterSize filterSigma filterWindow]
% Default: 'gaussian' [7 1 4]
% 
% TJ Young
% 20 February 2017

fsize = param(1); fsigma = param(2); fwindow = [param(3) param(3)];
filt = (fspecial(ftype,fsize,fsigma)); % Gaussian lowpass filter\
vv = medfilt2(vvIn,fwindow); % Median filtering in 2 directions
vvOut = conv2(vv,filt,'same'); % 2-D convolution with designed filte

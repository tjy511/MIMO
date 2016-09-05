function vvOut = peakConvol(vvIn, pk)

% vvOut = peakConvol(vvIn, pk)
%
% Peak Convolution Image Processing Method
%
% Input args:
% - vvIn = Pre-processed plane array (x,y) of power returns.
% - pk = Scalar value of +/- convolution.
%
% Output args: vvOut = Post-processed plane array (x,y) of power returns.
%
% Lai Bun Lok (edited by TJ Young)
% 08.12.2015

% Input parameters
vv = vvIn;

% ±4 peak convolution image processing technique
for n=1:length(vv);
    vvv=vv(:,n);
    [~,b]=max(vvv(5:length(vvv)-pk));
    b=b+pk;
    cf=vvv(b-pk:b+pk);
    www=conv(vvv,cf);
    yyy=spline(1:length(www),www,linspace(1,length(www),length(vvv)));
    vv(:,n)=yyy;
end
for n=1:length(vv)
    vvv=vv(n,:);
    [~,b]=max(vvv(5:length(vvv)-pk));
    b=b+pk;
    cf=vvv(b-pk:b+pk);
    www=conv(vvv,cf);
    yyy=spline(1:length(www),www,linspace(1,length(www),length(vvv)));
    vv(n,:)=yyy;
end

% Output parameters
vvOut = vv;
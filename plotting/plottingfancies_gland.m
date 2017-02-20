% Various ways to plot MIMO data
%
% TJ Young & Lai Bun Lok
% 08.12.2015

%% Load imagery file
load('array2d_20150703-1221.mat')
% array2d_20140506-1813.mat
% array2d_20140726-1727.mat
% array2d_20150703-1221.mat

%% Set depth parameters
depth = 635;%600;
depths = 630:650;%590:620;

% Correction to match Rs
depth = depth-10+1;
depths = depths-10+1; 

%% Plot 2D imagery

% X-axis
figure, hold on
surf(xxPix,Rs,20*log10(pp_slicex),'EdgeColor','none')
view(0,-90)
colormap(jet)
caxis([-100 -20])
legend = colorbar('Ticks',[-100 -80 -60 -40 -20]);
legend.Label.String = 'dB (Vrms)';
xlabel('Range (m)')
ylabel('Depth (m)')
title(['Vertical 2D profile (x-direction) at Date/Time: ', datestr(dateStamp)])

% Y-axis
figure,hold on
surf(yyPix,Rs,20*log10(pp_slicey),'EdgeColor','none')
view(0,-90)
colormap(jet)
caxis([-100 -20])
legend = colorbar('Ticks',[-100 -80 -60 -40 -20]);
legend.Label.String = 'dB (Vrms)';
xlabel('Range (m)')
ylabel('Depth (m)')
title(['Vertical 2D profile (y-direction) at Date/Time: ', datestr(dateStamp)])

%% Simple layer plot (pre-smooth)
figure, hold on
x = xxPix(depth,:);
y = yyPix(depth,:);
v = dB(imgPlane(:,:,depth));
surf(x,y,v,'EdgeColor','none')
colormap(jet)

%% Simple layer plot (post-smooth)
clear power 
% SR image post processing method
counter = 0; tt = 1;
for ii = depths
    counter = counter + 1;
    vv = peakConvol(imgPlane(:,:,ii),4);
    power(:,:,counter) = vv;
end

figure
surf(yPix,xPix,20*log10(abs(vv)))
title(['2D image at depth= ', num2str(depths),'metres', datestr(dateStamp(tt,1))])
colorbar
%xlim([-120 120])
%ylim([-120 120])
axis equal
%caxis([-100 -40])
%axis([-250 250 -250 250 -160 -40])
xlabel('x-position (m)')
ylabel('y-position (m)')
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure

%% Contour plot
clear v
% Store the matrices x, y, z, and v from the data set.
x = 1:100;
y = x;
z = flip(1:length(depths));
v = db(power);
% Limit to 9 planes in the x-z direction
Sx = 0;
Sy = 0;
Sz = 1:length(depths);

cvals = linspace(-60,-10,6); % Drawing contour lines

figure
contourslice(x,y,z,v,Sx,Sy,Sz,cvals)
xlabel('x'); ylabel('y');zlabel('z');
axis([0,100,0,100,z(end),z(1)])
colormap(jet)
legend = colorbar;
campos([0,-20,7])

%% Isosurface plot
thresh = -25; 
figure, hold on
p = patch(isosurface(x,y,z,v,thresh));
isonormals(x,y,z,v,p)
p.FaceColor = [0.5 0.5 0.5];
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight(250,-110)
lighting gouraud


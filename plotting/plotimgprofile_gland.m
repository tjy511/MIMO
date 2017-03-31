function [fig,h] = plotimgprofile_gland(xx,yy,zz)

% Quick function to plot 2-D depth profile of Greenland Store Data
%
% Inputs: 
% - Locations in x and y axes (normally 641x100 in size)
% - Magnitude (Amplitude) return (absolute) in z (depth) axis
%
% TJ Young
% 08.12.2015

%% Plot 2D imagery

% X-axis
fig = figure; 
hold on, box on, grid on
h = surf(xx,yy,db(zz,'voltage'),'EdgeColor','none');
view(0,-90)
colormap(jet)
caxis([-100 -20])
legend = colorbar('Ticks',[-100 -80 -60 -40 -20]);
legend.Label.String = 'dB (V_{RMS})';
xlabel('Cross-range [m]')
ylabel('Depth [m]')
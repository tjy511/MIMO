function plotimgdepth_gland(xx,yy,zz)

% Quick function to plot 2-D depth layer of Greenland Store Data
%
% Inputs: 
% - Locations in x and y axes (normally 100x100 in size)
% - Power return (absolute) in z (depth) axis
%
% TJ Young
% 08.12.2015

%% Plot 2D imagery

figure, hold on, axis equal
surf(xx,yy,zz,'EdgeColor','none')
view(0,90)
colormap(jet)
caxis([-100 -20])
legend = colorbar('Ticks',[-100 -80 -60 -40 -20]);
legend.Label.String = 'dB (Vrms)';
xlabel('X-position (m)')
ylabel('Y-position (m)')
%title(['2D image at depth= ', num2str(depth),'metres', datestr(dateStamp(tt,1))])
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
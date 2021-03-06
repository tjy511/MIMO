function fig = plotcumlayers_gland(x0,y0,z0,t,pxy)

% Function to plot layers in 3D
%
% [x,y,z] is a point on the layer.
%
% TJ Young
% 15 November 2016

% Establish normal vector of plane
%a = 0-x0;
%b = 0-y0;
%c = 0-z0;

% Plotting fancies
fig = figure; hold on, grid on
%plot3(0,0,0,'kp') % Location of radar unit
%quiver3(0,0,0,z0(end)*1.5,'Color',[0.6 0.6 0.6],'lineWidth',1.5) % Reference depth vector



%for cc = 1:size(x0,2)
    
    % Calculate plane
    x = pxy;
    [X,Y] = meshgrid(x);
    %Z = (a(cc).*X + b(cc).*Y) ./ -c(cc) + z0(cc);
    
    % Plot layer components iteratively
    %quiver3(0,0,0,x0(cc).*1.1,y0(cc).*1.1,z0(cc).*1.1,'k'); % Vector to identified peak
    
    plot3(x0,y0,t,'markerSize',10) % Location of identified peak
    %surf(X,Y,zeros(size(X)));
    if nargin == 5
        %z = squeeze(ppr(cc,:,:)); 
        %ax(cc) = patch(surf2patch(X,Y,Z,db(z)),'lineStyle','none'); % Plot layer with dB patched
    else
         % Plot layer with uniform patch as depth
    end
    
%end

% Figure graphics manipulator
%set(ax,'faceAlpha',0.5) % Transparancy
shading flat
legend = colorbar; 
%zlim([z0(end)-50 0])
view(0,90)
xlabel('x-position (m)'); 
ylabel('y-position (m)'); 
zlabel('Depth (m)');

if nargin == 5
    legend.Label.String = 'dB (Vrms)';
    colormap(jet)
    caxis([-100 -20])
else
    legend.Label.String = 'Depth (m)';
    caxis([z0(end) 0])
end
    
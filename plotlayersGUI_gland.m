%function fig = plotlayersGUI_gland(x0,y0,z0,pxy,ppr)
function PlotGUI(hObject,event)

%% Called by SimpleGUI to do the plotting
% hObject is the button and eventdata is unused.
%function PlotGUI(hObject,eventdata)
global Slider;% Slider is a handle to the slider.
global plane
global pxy

% Gets the value of the parameter from the slider.
Param = get(hObject,'Value');

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(Param),...
'Position', [460 55 60 20]);

x0 = plane.x0;
y0 = plane.y0;
z0 = plane.z0;

% Establish normal vector of plane
a = 0-x0;
b = 0-y0;
c = 0-z0;

    x = pxy(1,:);
    [X,Y] = meshgrid(x);
    Z = (a(1,1).*X + b(1,1).*Y) ./ -c(1,1) + z0(1,1);


hold on
plot3(0,0,0,'kp'); % Location of radar unit
quiver3(0,0,0,z0(1,end)*1.5,'Color',[0.6 0.6 0.6],'lineWidth',1.5) % Reference depth vector
ax = surf(X,Y,Z); % Plot layer with uniform patch as depth
hold off
set(ax,'xdata',X,'ydata',Y,'zdata',Z);

% %Plots the Graph.
% x=linspace(0,10,1000);
% k = Param;
% y = sin(k*x);
% %h1 = plot(x,zeros(1,1000));
% 
% h1 = plot(0,0); % dummy vector
% hold on
% h2 = plot(x,y);
% h3 = plot(x,y);
% h4 = plot(x,y);
% h5 = plot(x,y);
% hold off
% set(h2,'xdata',x,'ydata',2*y);
% set(h3,'xdata',x,'ydata',3*y);
% set(h4,'xdata',x,'ydata',4*y);
% set(h5,'xdata',x,'ydata',5*y);

%{
% Function to plot layers in 3D
%
% [x,y,z] is a point on the layer.
%
% TJ Young
% 15 November 2016

% Establish normal vector of plane
a = 0-x0;
b = 0-y0;
c = 0-z0;

% Plotting fancies
fig = figure; hold on, grid on
plot3(0,0,0,'kp') % Location of radar unit
quiver3(0,0,0,z0(end)*1.5,'Color',[0.6 0.6 0.6],'lineWidth',1.5) % Reference depth vector

for cc = 1:size(x0,2)
    
    % Calculate plane
    x = pxy(cc,:);
    [X,Y] = meshgrid(x);
    Z = (a(cc).*X + b(cc).*Y) ./ -c(cc) + z0(cc);
    
    % Plot layer components iteratively
    quiver3(0,0,0,x0(cc).*1.1,y0(cc).*1.1,z0(cc).*1.1,'k'); % Vector to identified peak
    plot3(x0(cc),y0(cc),z0(cc),'k*') % Location of identified peak
    if nargin == 5
        z = squeeze(ppr(cc,:,:)); 
        ax(cc) = patch(surf2patch(X,Y,Z,db(z)),'lineStyle','none'); % Plot layer with dB patched
    else
        ax(cc) = surf(X,Y,Z); % Plot layer with uniform patch as depth
    end
    
end

% Figure graphics manipulator
set(ax,'faceAlpha',0.5) % Transparancy
shading flat
legend = colorbar; 
zlim([z0(end)-50 0])
view(-15,15)
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
%}
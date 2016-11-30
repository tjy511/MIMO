%% Called by SimpeGUI to do the plotting
% hObject is the button and eventdata is unused.
function PlotGUI(hObject,eventdata)
global Slider;% Slider is a handle to the slider.

% Gets the value of the parameter from the slider.
Param = get(hObject,'Value');

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(Param),...
'Position', [460 55 60 20]);

% Plots the Graph.
x=linspace(0,10,1000);
k = Param;
y = sin(k*x);
%h1 = plot(x,zeros(1,1000));

h1 = plot(0,0); % dummy vector
hold on
h2 = plot(x,y);
h3 = plot(x,y);
h4 = plot(x,y);
h5 = plot(x,y);
hold off
set(h2,'xdata',x,'ydata',2*y);
set(h3,'xdata',x,'ydata',3*y);
set(h4,'xdata',x,'ydata',4*y);
set(h5,'xdata',x,'ydata',5*y);
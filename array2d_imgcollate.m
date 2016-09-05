% Script to make movies out of 2D processing data

%% Parameters

% Get list of processed files
fileloc = '/Users/tjy22/Documents/School/PhD/radar/data/process/mimo/unattended/';
data_list=dir(fullfile(fileloc,'*.mat'));

% X or Y
xy = 'yyy'; % xxx or yyy

%% Create M

hfig = figure;

%for frame = 1
for frame = 1:length(data_list)
    
    % Load data file
    load(data_list(frame).name)
    
    % Switch between x and y axes
    switch xy
        case 'xxx'
            xData = xxPix;
            yData = Rs;
            zData = pp_slicex;
            stringData = 'x';
        case 'yyy'
            xData = yyPix;
            yData = Rs;
            zData = pp_slicey;
            stringData = 'y';
    end
    
    % Plot graph
    surf(xData,yData,20*log10(zData),'EdgeColor','none')
    % Set axes
    view(0,-90)
    %axis([-150 150 0 300])
    % Set legend
    colormap(jet)
    caxis([-100 -20])
    legend = colorbar('Ticks',[-100 -80 -60 -40 -20]);
    legend.Label.String = 'dB (Vrms)';
    % Label graph
    xlabel('Range (m)')
    ylabel('Depth (m)')
    title(['Vertical 2D profile (', stringData', '-direction)'])
    S1 = strcat(['Date/Time: ', datestr(dateStamp)]);
    text(50,25,S1)
    %text(25,10,S1)
    
    % Assign graph to array
    M(frame) = getframe(hfig);
end


%% Save frames
switch xy
    case 'xxx'
        Mx = M;
    case 'yyy'
        My = M;
end

%% Play Movie

FigHandle = figure;
axis off
set(FigHandle, 'Position', [100, 100, 650, 500]);
movie(My,1,3);

%% Save movie
switch xy
    case 'xxx'
        label = 'Mx';
    case 'yyy'
        label = 'My';
end
videoName = strcat(label,'.mp4');
v = VideoWriter(videoName, 'MPEG-4'); % Create video handle
v.FrameRate = 3; 
open(v) % Open video handle
writeVideo(v, My); % Save v to handle; CHANGE VARIABLE AS NEEDED
disp(['Saving movie as file: ', videoName])
close(v) % Close video handle
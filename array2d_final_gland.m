% Script to display the virtual element positions in 2D mimo array
% Modified for L shape
% Script written by Lai Bun Lok


%% Housekeeping

clear all
close all

plotelements=0; % Plot (1) out elements or not (0)

%% Antenna separation

%dphy=0.72 % Physical antenna separation (in metres), centre-to-centre
dphy=0.83; % Greenland setting

% Position of Tx element real positions (m)

%xt1=2; % offset from origin
xt1=2.49; % Greenland setup
yt1=0;

xt2=xt1+dphy;
yt2=0;

xt3=xt2+dphy;
yt3=0;

xt4=xt3+dphy;
yt4=0;

xt5=xt4+dphy;
yt5=0;

xt6=xt5+dphy;
yt6=0;

xt7=xt6+dphy;
yt7=0;

xt8=xt7+dphy;
yt8=0;

%xt9=xt8+dphy;
%yt9=0;

tx=zeros(8,2);
tx(1,1)=xt1; tx(1,2)=yt1;
tx(2,1)=xt2; tx(2,2)=yt2;
tx(3,1)=xt3; tx(3,2)=yt3;
tx(4,1)=xt4; tx(4,2)=yt4;
tx(5,1)=xt5; tx(5,2)=yt5;
tx(6,1)=xt6; tx(6,2)=yt6;
tx(7,1)=xt7; tx(7,2)=yt7;
tx(8,1)=xt8; tx(8,2)=yt8;
%tx(9,1)=xt9; tx(9,2)=yt9;

% Position of Rx element real positions (m)

xr1=0;
%yr1=-2; % offset from origin
yr1=-2.49; % Greenland set up

xr2=xr1;
yr2=yr1-dphy;

xr3=xr1;
yr3=yr2-dphy;

xr4=xr1;
yr4=yr3-dphy;

xr5=xr1;
yr5=yr4-dphy;

xr6=xr1;
yr6=yr5-dphy;

xr7=xr1;
yr7=yr6-dphy;

xr8=xr1;
yr8=yr7-dphy;

%xr9=xr1;
%yr9=yr8-dphy;

rx=zeros(8,2);
rx(1,1)=xr1; rx(1,2)=yr1;
rx(2,1)=xr2; rx(2,2)=yr2;
rx(3,1)=xr3; rx(3,2)=yr3;
rx(4,1)=xr4; rx(4,2)=yr4;
rx(5,1)=xr5; rx(5,2)=yr5;
rx(6,1)=xr6; rx(6,2)=yr6;
rx(7,1)=xr7; rx(7,2)=yr7;
rx(8,1)=xr8; rx(8,2)=yr8;
%rx(9,1)=xr9; rx(9,2)=yr9;

% Calculate the virtual element positions

for x=1:8
    ve(x,1)=tx(1,1)+(rx(x,1)-tx(1,1))/2;
    ve(x,2)=rx(x,2)+(tx(1,2)-rx(x,2))/2;
end


for x=9:16
    ve(x,1)=tx(2,1)+(rx(x-8,1)-tx(2,1))/2;
    ve(x,2)=rx(x-8,2)+(tx(2,2)-rx(x-8,2))/2;
end


for x=17:24
    ve(x,1)=tx(3,1)+(rx(x-16,1)-tx(3,1))/2;
    ve(x,2)=rx(x-16,2)+(tx(3,2)-rx(x-16,2))/2;
end

for x=25:32
    ve(x,1)=tx(4,1)+(rx(x-24,1)-tx(4,1))/2;
    ve(x,2)=rx(x-24,2)+(tx(4,2)-rx(x-24,2))/2;
end

for x=33:40
    ve(x,1)=tx(5,1)+(rx(x-32,1)-tx(5,1))/2;
    ve(x,2)=rx(x-32,2)+(tx(5,2)-rx(x-32,2))/2;
end

for x=41:48
    ve(x,1)=tx(6,1)+(rx(x-40,1)-tx(6,1))/2;
    ve(x,2)=rx(x-40,2)+(tx(6,2)-rx(x-40,2))/2;
end

for x=49:56
    ve(x,1)=tx(7,1)+(rx(x-48,1)-tx(7,1))/2;
    ve(x,2)=rx(x-48,2)+(tx(7,2)-rx(x-48,2))/2;
end

for x=57:64
    ve(x,1)=tx(8,1)+(rx(x-56,1)-tx(8,1))/2;
    ve(x,2)=rx(x-56,2)+(tx(8,2)-rx(x-56,2))/2;
end

% for x=73:81
%     ve(x,1)=tx(9,1)+(rx(x-72,1)-tx(9,1))/2;
%     ve(x,2)=rx(x-72,2)+(tx(9,2)-rx(x-72,2))/2;
% end

%% Plot them out in xy

if plotelements==1
    figure
    hold
    xlabel('x position (m)')
    ylabel('y position (m)')
    
    % Plot Tx and Rx antennas
    for k=1:8
        plot(tx(k,1),tx(k,2),'+')
        plot(rx(k,1),rx(k,2),'o')
        text(tx(k,1)-0.08,tx(k,2)+0.4,num2str(k))
        text(rx(k,1)-0.4,rx(k,2)-0.04,num2str(k))
    end
    text(tx(1,1)-dphy-0.08,tx(1,2)+0.4,'Tx')
    text(rx(1,1)-0.4,rx(1,2)+dphy+0.04,'Rx')
    
    % Restrict axes
    xlim([-1 10])
    ylim([-10 1])
    
    % Plot virtual antennas
    for k=1:64
        plot(ve(k,1),ve(k,2),'rx')
    end
    
end
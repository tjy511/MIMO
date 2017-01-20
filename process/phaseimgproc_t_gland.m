% Function to process MIMO data into 2-dimensional slices through depth. 
%
% Outputs:
% - xxPix/yyPix: Matrix of locations of cone in 2D space (x or y)
% - Rs: Matrix lof locations on cone in 3rd space (z)
% - pp_slicex/y: Power return at x or y direction
% - imgPlane: Power return in all 3 dimensions
%
% Written by Lai Bun Lok
% Edited by TJ Young
% Last edited: 05/09/2016

%% Prep work and identify parameters

% Identify data and type
startup
deployment = 3; % 1 2 3
processing = 'single'; % Single or multiple chirps 'ts' 'single'
%fig_dir = strcat(rwd,'/data/process/mimo/attended/');
fig_dir = '~/Downloads/';

bstart = 1; % Starting burst (normally 1)
leapFrog = 1*24; % In bursts. Set to 1 if processing = 'single'.

% Parameters
dphy = 0.83; % Antenna separation (centre-to-centre)
fs=40000; % Samples per chirp (Note 40001 for deployment4)
Npix=100; % number of pixels in Npix*Npix image plane.
fov=25; % field of view plus-minus degrees.
%R=617.23; % Bed depth (Greenland)
Rs=10:1:650; % Range slices

% Pre-allocate arrays
imgPlane=zeros(Npix,Npix,length(Rs));
xxPix=zeros(length(Rs),Npix);
yyPix=zeros(length(Rs),Npix);
pp_slicey=zeros(length(Rs),Npix);
pp_slicex=zeros(length(Rs),Npix);
r=zeros(Npix*Npix,64);

%% Load antenna locations
folderLoc = '/Users/tjy511/OneDrive - University Of Cambridge/radar/data/field/array/';
switch processing
    case 'ts'
        switch deployment
            case 1
                [txLoc,rxLoc,ve] = antennaLoc('store1',1);
                tsList = readtable('radarlist1.dat');
                data_dir = strcat(folderLoc,'unattended/deployment1/timeseries/');
            case 2
                [txLoc,rxLoc,ve] = antennaLoc('store2',1);
                tsList = readtable('radarlist3a.dat');
                data_dir = strcat(folderLoc,'unattended/deployment3/timeseries/');
            case 3
                [txLoc,rxLoc,ve] = antennaLoc('store3',1);
                tsList = readtable('radarlist4a.dat');
                data_dir = strcat(folderLoc,'unattended/deployment4/timeseries/');
        end
    case 'single'
        switch deployment
            case 1
                [txLoc,rxLoc,ve] = antennaLoc('store1',1);
                tsList = 'Survey_2014-05-06_181256.dat';
                data_dir = strcat(folderLoc,'attended/combined/20140506/');
            case 2
                [txLoc,rxLoc,ve] = antennaLoc('store2',1);
                tsList = 'Survey-2014-07-26_172711.dat';
                data_dir = strcat(folderLoc,'attended/combined/20140726/');
            case 3
                [txLoc,rxLoc,ve] = antennaLoc('store3',1);
                tsList = 'Survey_2015-07-03_122127.dat';
                data_dir = strcat(folderLoc,'attended/combined/20150705/');
        end
end
ant = [txLoc(:,1) rxLoc(:,2)];

%% Load file

for fileNum = bstart:leapFrog:size(tsList,1)
    
    % Read in the raw data
    switch processing
        case 'ts'
            xf = fileNum;
            bnum = tsList.burst(xf);
            vdat = LoadBurstRMB4(fullfile(data_dir,tsList.file{xf}),bnum,fs);
        case 'single'
            vdat = LoadBurstRMB4(fullfile(data_dir,tsList),1,fs);
    end
    dateStamp=vdat.TimeStamp;
    fileDate = datestr(dateStamp,'yyyymmdd-HHMM');
    fileName = strcat('new_array2d_', fileDate, '.mat');
    
    if exist(fileName,'file')==0
    disp('*********************************************************')
    disp(['Running script on burst obtained at: ',datestr(dateStamp)])
    
    %% Iterate processing for every depth step
    for ss=1:length(Rs)
        
        R=Rs(ss); % Range at step ss
        ImgPlane=zeros(Npix*Npix,2); % Pre-allocate image plane
        disp(['Processing depth step at: ', num2str(R), ' m'])
        
        % Pixel size across depth image plane
        dx=R*tan(fov*pi/180)*2/Npix;
        dy=dx;
        
        % Create a 1xNpix matrix containing locations of each pixel on the image plane
        for x=1:Npix
            for y=1:Npix
                ImgPlane(y+(x-1)*Npix,1)=(-Npix/2)*dx+(x-1)*dx+dx/2;
                ImgPlane(y+(x-1)*Npix,2)=(-Npix/2)*dy+(y-1)*dy+dy/2;
            end
        end
        
%         % Plot image slices per depth
%         figure
%         hold all
%         for k=1:Npix*Npix
%             plot(ImgPlane(k,1),ImgPlane(k,2), 'b+')
%         end
%         xlabel('x position (m)')
%         ylabel('y position (m)')
%         title(['Pixel locations at range R=', num2str(R), ' metres'])
%         Recentre virtual element array to origin
        ve_offset=ve-(1+7*dphy/4);
        ve_offset(:,2)=-1*(ve_offset(:,2)+2*(1+7*dphy/4));
%         plot(ve_offset(:,1),ve_offset(:,2),'d')
        
        % Calculate distance r of virtual elements to each pixel on image plane
        for k=1:Npix*Npix
            for v=1:length(ve)
                r(k,v)=sqrt((ImgPlane(k,1)-ve_offset(v,1))^2+(ImgPlane(k,2)-ve_offset(v,2))^2+R^2);
            end
        end
        
        %% Work out the total date files to process
        
        % Apply the calibrated range processing algorithms
        range64_calib_gland % Correct for cable length
        lambda_g=3e8/sqrt(er)/fc; %
        
        P_pix=zeros(Npix*Npix,1); % Pre-allocate array
        
        for k=1:Npix*Npix
            for v=1:length(ve)
                
                % Select the range bin corresponding to r for each virtual array element
                B(k,v)=round(r(k,v)*padfactor/deltaR);
                
                % Weight value of the return in this range bin
                % Not sure what "weight" means tbh.
                weight=exp(-j*4*pi*r(k,v)/lambda_g);
                specWeighted=ampspec(v,B(k,v))*weight;
                
                % Sum for all virtual elements
                % Not sure what this is again...
                %         P_pix(k,1)=P_pix(k,1)+ampspec(v,B(k,v));
                P_pix(k,1)=P_pix(k,1)+specWeighted;
            end
        end
        
        % Reorganise pixel matrix into an Npix*Npix array
        PP_pix=reshape(P_pix,Npix,Npix);

        %         % Visualise results
        close all
        %         %figure
        yPix=ImgPlane(1:Npix,2);
        xPix=yPix;
        %
        %         surf(yPix,xPix,20*log10(abs(PP_pix)))
        %         %view([0 -90])
        %         title(['2D image at range R=', num2str(R),'metres', datestr(dateStamp(tt,1))])
        %         colorbar
        %         %xlim([-120 120])
        %         %ylim([-120 120])
        %         axis equal
        %         xlabel('x-position (m)')
        %         ylabel('y-position (m)')
        %         %caxis([-40 -5])
        %         %caxis([-30 -10])
        %         %caxis([0 0.5]) % Site 1 linear scale
        %         %caxis([0 0.2]) % Site 3 linear scale
        %         set(gcf,'Renderer','Zbuffer')
        %         figure
        %         % SR image post processing method
        %         % % vv=abs(PP_pix);
        %         % % for n=1:length(vv);
        %         % %     vvv=vv(:,n);
        %         % %     [a,b]=max(vvv(5:length(vvv)-4));
        %         % %     b=b+4;
        %         % %     cf=vvv(b-4:b+4);
        %         % %     www=conv(vvv,cf);
        %         % %     yyy=spline(1:length(www),www,linspace(1,length(www),length(vvv)));
        %         % %     vv(:,n)=yyy;
        %         % % end
        %         % %
        %         % % for n=1:length(vv)
        %         % %     vvv=vv(n,:);
        %         % %     [a,b]=max(vvv(5:length(vvv)-4));
        %         % %     b=b+4;
        %         % %     cf=vvv(b-4:b+4);
        %         % %     www=conv(vvv,cf);
        %         % %     yyy=spline(1:length(www),www,linspace(1,length(www),length(vvv)));
        %         % %     vv(n,:)=yyy;
        %         % % end
        %         % % surf(yPix,xPix,20*log10(abs(vv)))
        %         % % title(['2D image at range R=', num2str(R),'metres', datestr(dateStamp(tt,1))])
        %         % % colorbar
        %         % % %xlim([-120 120])
        %         % % %ylim([-120 120])
        %         % % axis equal
        %         % % %caxis([-100 -40])
        %         % % %axis([-250 250 -250 250 -160 -40])
        %         % % xlabel('x-position (m)')
        %         % % ylabel('y-position (m)')
        %         % % set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        %
        %         %M(idxIm)=getframe(gcf);
        %         %idxIm=idxIm+1;
        %         %close gcf
        % %    end
        
        imgPlane(:,:,ss)=abs(PP_pix);
        yyPix(ss,:)=yPix;
        xxPix(ss,:)=xPix;
        pp_slicex(ss,:)=(abs(PP_pix(Npix/2,:)));
        pp_slicey(ss,:)=(abs(PP_pix(:,Npix/2)));
    end
    
    % Save file
    fileDate = datestr(dateStamp,'yyyymmdd-HHMM');
    fileName = strcat('array2d_', fileDate, '.mat');
    fileNameFull = strcat(fig_dir,fileName);
    save(fileNameFull);
    disp(['Saved data results as file: ', fileName])
    end
end
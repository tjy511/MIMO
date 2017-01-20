% Quick script to automate profile recreation through different antenna
% combinations

%%  Set 1.1

txrx = [32 1]; dPhy = [0.25 0.5 1 2]; 
ant = {'isotropic' 'pencil' 'bowtie' 'helix' 'helix2' 'helix3' 'helix4' 'helix5' 'dipole'};
for ii = 1:length(dPhy)
    for jj = 1:length(ant)
        disp(['Running script for parameters: ',ant{jj},' || ',num2str(txrx),' || ',num2str(dPhy(ii))])
        profile_recreation(ant{jj},txrx,dPhy(ii))
        disp(' ')
        close all
    end
end

%% Set 1.2
txrx = [8 8]; dPhy = [0.5 1 1.5 2]; 
ant = {'isotropic' 'pencil' 'bowtie' 'helix' 'helix2' 'helix3' 'helix4' 'helix5' 'dipole'};
for ii = 1:length(dPhy)
    for jj = 1:length(ant)
        disp(['Running script for parameters: ',ant{jj},' || ',num2str(txrx),' || ',num2str(dPhy(ii))])
        profile_recreation(ant{jj},txrx,dPhy(ii))
        disp(' ')
        close all
    end
end

%% Set 1.2.0
txrx = [24 1]; dPhy = [0.5 1 1.5 2]; 
ant = {'isotropic' 'pencil' 'bowtie' 'helix' 'dipole'}; 
for ii = 1:length(dPhy)
    for jj = 1:length(ant)
        disp(['Running script for parameters: ',ant{jj},' || ',num2str(txrx),' || ',num2str(dPhy(ii))])
        profile_recreation(ant{jj},txrx,dPhy(ii))
        disp(' ')
        close all
    end
end

%% Set 1.3
txrx = [8 1]; dPhy = [0.5 1 1.5 2]; 
ant = {'isotropic' 'pencil' 'bowtie' 'helix' 'dipole'};
for ii = 1:length(dPhy)
    for jj = 1:length(ant)
        disp(['Running script for parameters: ',ant{jj},' || ',num2str(txrx),' || ',num2str(dPhy(ii))])
        profile_recreation(ant{jj},txrx,dPhy(ii))
        disp(' ')
        close all
    end
end

%% Set 2.1
txrx = [8 8]; dPhy = [0.5 1 1.5 2]; 
ant = {'isotropic' 'pencil' 'bowtie' 'helix' 'dipole'};
for ii = 1:length(dPhy)
    for jj = 1:length(ant)
        disp(['Running script for parameters: ',ant{jj},' || ',num2str(txrx),' || ',num2str(dPhy(ii))])
        profile_recreation(ant{jj},txrx,dPhy(ii))
        disp(' ')
        close all
    end
end

%% Set 2.2
txrx = [12 4]; dPhy = [0.5 1 1.5 2]; 
ant = {'isotropic' 'pencil' 'bowtie' 'helix' 'dipole'};
for ii = 1:length(dPhy)
    for jj = 1:length(ant)
        disp(['Running script for parameters: ',ant{jj},' || ',num2str(txrx),' || ',num2str(dPhy(ii))])
        profile_recreation(ant{jj},txrx,dPhy(ii))
        disp(' ')
        close all
    end
end

%% Set 2.3
txrx = [14 2]; dPhy = [0.5 1 1.5 2]; 
ant = {'isotropic' 'pencil' 'bowtie' 'helix' 'dipole'};
for ii = 1:length(dPhy)
    for jj = 1:length(ant)
        disp(['Running script for parameters: ',ant{jj},' || ',num2str(txrx),' || ',num2str(dPhy(ii))])
        profile_recreation(ant{jj},txrx,dPhy(ii))
        disp(' ')
        close all
    end
end
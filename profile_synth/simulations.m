[theta,RA,RE] = antennaBP('pencil',1,0.00333);
hp_bw = hpbw(theta,RE);

% a = 0.00148675 == HPBW of 45.001
% a = 0.00333 == HPBW of 30.0680
% a = 0.0075 == HPBW of 20.0360
% a = 0.0301 == HPBW of 10.0020
% a = 0.1205 == HPBW of 5.000

%% Experiment 1

[fig1,h1] = profileRecreation(mimo_config_sim11);
[fig2,h2] = profileRecreation(mimo_config_sim12);
[fig3,h3] = profileRecreation(mimo_config_sim13);
[fig4,h4] = profileRecreation(mimo_config_sim14);
%% Experiment 2

[fig1,h1] = profileRecreation(mimo_config_sim21);
[fig2,h2] = profileRecreation(mimo_config_sim22);
[fig3,h3] = profileRecreation(mimo_config_sim23);
[fig4,h4] = profileRecreation(mimo_config_sim24);

%% Experiment 3

[fig1,h1] = profileRecreation(mimo_config_sim31);
[fig2,h2] = profileRecreation(mimo_config_sim32);
[fig3,h3] = profileRecreation(mimo_config_sim33);
[fig4,h4] = profileRecreation(mimo_config_sim34);

%% Experiment 4

[fig1,h1] = profileRecreation(mimo_config_sim41);
[fig2,h2] = profileRecreation(mimo_config_sim42);
[fig3,h3] = profileRecreation(mimo_config_sim43);
[fig4,h4] = profileRecreation(mimo_config_sim44);

%% Simulation of radar setup

% Make sure antLoc in script is set to store1
profile_recreation('bowtie',[8 8],1,1); % Last 2 terms are dummy
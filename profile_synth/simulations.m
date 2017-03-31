%[theta,RA,RE] = antennaBP('pencil',1,0.00333);
%hp_bw = hpbw(theta,RE);

% a = 0.00148675 == HPBW of 45.001
% a = 0.00333 == HPBW of 30.0680
% a = 0.0075 == HPBW of 20.0360
% a = 0.0301 == HPBW of 10.0020
% a = 0.1205 == HPBW of 5.000

%% Field simulations
[fig1,h1] = profileRecreation(mimo_config_field_mtn);

%% Simulation 1

[fig11,h11] = profileRecreation(mimo_config_sim11);
[fig12,h12] = profileRecreation(mimo_config_sim12);
[fig13,h13] = profileRecreation(mimo_config_sim13);
[fig14,h14] = profileRecreation(mimo_config_sim14);
%% Simulation 2

[fig21,h21] = profileRecreation(mimo_config_sim21);
[fig22,h22] = profileRecreation(mimo_config_sim22);
[fig23,h23] = profileRecreation(mimo_config_sim23);
[fig24,h24] = profileRecreation(mimo_config_sim24);

%% Simulation 3

[fig31,h31] = profileRecreation(mimo_config_sim31);
[fig32,h32] = profileRecreation(mimo_config_sim32);
[fig33,h33] = profileRecreation(mimo_config_sim33);
[fig34,h34] = profileRecreation(mimo_config_sim34);

%% Simulation of radar setup

% Make sure antLoc in script is set to store1
%profile_recreation('bowtie',[8 8],1,1); % Last 2 terms are dummy

%% Saving
cd('/Users/tjy511/Google Drive/Academic/papers/paper3/figs/simulation/field_setup/');
export_fig(fig1,'bowtie_8-8_0.83_PR_2d.png','-m6');
%%
cd('/Users/tjy511/Google Drive/Academic/papers/paper3/figs/simulation/experiments/experiment1');
export_fig(fig11,'pencil_10_1_32-1_PR_2d.png','-m6');
export_fig(fig12,'pencil_20_1_32-1_PR_2d.png','-m6');
export_fig(fig13,'pencil_30_1_32-1_PR_2d.png','-m6');
export_fig(fig14,'pencil_45_1_32-1_PR_2d.png','-m6');
%%
cd('/Users/tjy511/Google Drive/Academic/papers/paper3/figs/simulation/experiments/experiment2');
export_fig(fig21,'pencil_30_0.25_32-1_PR_2d.png','-m6');
export_fig(fig22,'pencil_30_0.5_32-1_PR_2d.png','-m6');
export_fig(fig23,'pencil_30_1_32-1_PR_2d.png','-m6');
export_fig(fig24,'pencil_30_2_32-1_PR_2d.png','-m6');
%%
cd('/Users/tjy511/Google Drive/Academic/papers/paper3/figs/simulation/experiments/experiment3');
export_fig(fig31,'pencil_30_1_32-1_PR_2d.png','-m6');
export_fig(fig32,'pencil_30_1_24-1_PR_2d.png','-m6');
export_fig(fig33,'pencil_30_1_16-1_PR_2d.png','-m6');
export_fig(fig34,'pencil_30_1_8-1_PR_2d.png','-m6');
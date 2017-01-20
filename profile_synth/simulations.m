[theta,RA,RE] = antennaBP('pencil',1,0.00333);
hp_bw = hpbw(theta,RE);

% a = 0.00148675 == HPBW of 45.001
% a = 0.00333 == HPBW of 30.0680
% a = 0.0075 == HPBW of 20.0360
% a = 0.0301 == HPBW of 10.0020
% a = 0.1205 == HPBW of 5.000

%% Experiment 1

profile_recreation('pencil',[32 1],1,0.1205)
profile_recreation('pencil',[32 1],1,0.0301)
profile_recreation('pencil',[32 1],1,0.0075)
profile_recreation('pencil',[32 1],1,0.00333)
profile_recreation('pencil',[32 1],1,0.00148675)

%% Experiment 2
profile_recreation('pencil',[32 1],1/4,0.00333)
profile_recreation('pencil',[32 1],1/2,0.00333)
profile_recreation('pencil',[32 1],1,0.00333)
profile_recreation('pencil',[32 1],2,0.00333)


%% Experiment 3

profile_recreation('pencil',[32 1],1,0.00333)
profile_recreation('pencil',[24 1],1,0.00333)
profile_recreation('pencil',[16 1],1,0.00333)
profile_recreation('pencil',[8 1],1,0.00333)
profile_recreation('pencil',[8 8],1,0.00333)

%% Experiment 4

profile_recreation('isotropic',[32 1],1,0.00333)
profile_recreation('pencil',[32 1],1,0.00333)
profile_recreation('bowtie',[32 1],1,0.00333)
profile_recreation('helix',[32 1],1,0.00333)

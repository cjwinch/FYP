addpath('/Users/chris/Documents/FYP/MATLAB/cjw_pls_func');
xvars = {'meanRR', 'stdRR',    ...
       'meanHR', 'stdHR',    ...
       'RMSSD', 'SDSD',      ...
       'NN50', 'PNN50',      ...
       'VLF', 'HF', 'LF',    ...
       'PTOT', 'LFN', 'HFN', ...
       'LFHF',               ...
       'SampEn1', 'SampEn2', ...
       'meanIBI', 'stdRR', ...
       'meanResp', 'stdResp'
       };

load('participant1results.mat');
p1 = all_stats;
load('participant2results.mat');
p2 = all_stats;
load('participant3results.mat');
p3 = all_stats;
load('participant4results.mat');
p4 = all_stats;
load('participant6results.mat');
p6 = all_stats;
load('participant7results.mat');
p7 = all_stats;
load('participant8results.mat');
p8 = all_stats;
load('participant9results.mat');
p9 = all_stats;
   
m1 = p1(1,1);
m2 = p2(1,1);
m3 = p3(1,1);
m4 = p4(1,1);
m6 = p6(1,1);
m7 = p7(1,1);
m8 = p8(1,1);
m9 = p9(1,1);

%% Try to normalise the variables 



%    
%% Regression on participant 9 against pereceived stress

load('participant9results.mat');
use = [(1:8) (10:21)]; % Exclude VLF
X = all_stats(2:5,use);

Y = [6; 8; 9; 9];

model = cjw_pls2(X,Y,4);

% Biplot for X
p = [model.p];
t = [model.t];
biplot(p(:,1:2), 'scores', t(:,1:2), ...
       'varlabels', xvars(use), 'ObsLabels', {'Control', '1', '2', '3'});

% Biplot for Y
figure(2)
q = [model.q];
u = [model.u];
biplot(q(:,1:2), 'scores', u(:,1:2), ...
       'varlabels', {'Stress'}, 'ObsLabels', {'Control', '1', '2', '3'});

% Inner relation between X and Y
figure(3);
scatter(t(:,1), u(:,1), 'r');

% Scree plot
figure(4);
plot(cumsum([model.R2X])); hold on;
plot(cumsum([model.R2Y]));

% Try with PCA
[coeff,score,latent,tsquared] = pca(X);
biplot(coeff(:,1:2), 'scores', score(:,1:2), ...
       'varlabels', xvars(use), 'ObsLabels', {'Control', '1', '2', '3'});
   
%% Regression on participant 9 against pereceived stress and difficulty

% load('participant9results.mat');
X = all_stats(2:5,:);

Y = [[6; 8; 9; 9] [3; 8; 8; 9]];

model = cjw_pls2(X,Y,4);

% Biplot for X
p = [model.p];
t = [model.t];
biplot(p(:,1:2), 'scores', t(:,1:2), ...
       'varlabels', xvars, 'ObsLabels', {'Control', '1', '2', '3'});

% Biplot for Y
figure(2)
q = [model.q];
u = [model.u];
biplot(q(:,1:2), 'scores', u(:,1:2), ...
       'varlabels', {'Stress', 'Difficulty'}, 'ObsLabels', {'Control', '1', '2', '3'});
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

%    %% Abdi
%    X = all_stats(2:5,:);
% 
% Y = [6; 8; 9; 9];
% 
% 
%    
%% Regression on participant 9 against pereceived stress

load('participant9results.mat');

use = [(1:8) (10:18)]; % Exclude VLF, stdResp, second stdRR, meanResp

X = all_stats(1:5,use);
xvu = xvars(use);
Y = [4; 6; 8; 9; 9];

model = cjw_pls2(X,Y,5);

test =  [1 X(1,:)] * model(1).b_coeff_star; % should equal 4

% Biplot for X
p = [model.p];
t = [model.t];
biplot(p(:,1:2), 'scores', t(:,1:2), ...
       'varlabels', xvu, 'ObsLabels', {'Rest', 'Control', '1', '2', '3'});

xlabel('Latent Variable 1');
xlabel('Latent Variable 2');
xlim(2.6*[-1 1]);
ylim(2.6*[-1 1]);
pause
% CJWPlotV4('MIST_Participant9_X_Biplot_DNM' );  
   
% Biplot for Y
figure(2)
q = [model.q];
u = [model.u];
biplot(q(:,1:2), 'scores', u(:,1:2), ...
       'varlabels', {'Stress'}, 'ObsLabels', {'Rest', 'Control', '1', '2', '3'});

% Inner relation between X and Y
figure(3);
scatter(t(:,1), u(:,1), 'r');

% Scree plot
figure(4);
plot(100*cumsum([model.R2X])); hold on;
plot(100*cumsum([model.R2Y]));
ylabel('Variance Explained (\%)');
xlabel('Number of Components');
legend('$R^2_X$','$R^2_Y$', 'location', 'southeast');  
xlim([1 5]);
set(gca, 'XTick', 1:5);
ylim([0 100]);
CJWPlotV4('MIST_Participant9_X_Scree');

% Try with PCA
[coeff,score,latent,tsquared] = pca(zscore(X));
biplot(coeff(:,1:2), 'scores', score(:,1:2), ...
       'varlabels', xvars(use), 'ObsLabels', {'Control', '1', '2', '3'});
   
%% Regression on participant 9 against pereceived stress and difficulty

% load('participant9results.mat');
% X = all_stats(2:5,:);
% 
% Y = [[6; 8; 9; 9] [3; 8; 8; 9]];
% 
% model = cjw_pls2(X,Y,4);
% 
% % Biplot for X
% p = [model.p];
% t = [model.t];
% biplot(p(:,1:2), 'scores', t(:,1:2), ...
%        'varlabels', xvars, 'ObsLabels', {'Control', '1', '2', '3'});
% 
% % Biplot for Y
% figure(2)
% q = [model.q];
% u = [model.u];
% biplot(q(:,1:2), 'scores', u(:,1:2), ...
%        'varlabels', {'Stress', 'Difficulty'}, 'ObsLabels', {'Control', '1', '2', '3'});
%   
%% VIP
vipp = cjw_vip( X, Y, model );
% VIP=vip(X,Y,[model.t],[model.w],[model.q])
reduction_criteria = vipp>1 | (vipp<1 & abs(model(1).b_coeff)'>0.1);
reduced_model = cjw_pls2(X(:,reduction_criteria),Y,5);

p_r = [reduced_model.p];
t_r = [reduced_model.t];
biplot(p_r(:,1:2), 'scores', t_r(:,1:2), ...
       'varlabels', xvu(reduction_criteria), 'ObsLabels', {'Rest', 'Control', '1', '2', '3'});
xlabel('Latent Variable 1');
xlabel('Latent Variable 2');
xlim(2.6*[-1 1]);
ylim(2.6*[-1 1]);
pause
CJWPlotV4('MIST_Participant9_X_Biplot_Reduced_Model_DNM' );  
   
reduction_criteria2 = vipp>1 | (abs(model(1).b_coeff)'>0.1);
reduced_model2 = cjw_pls2(X(:,reduction_criteria2),Y,5);



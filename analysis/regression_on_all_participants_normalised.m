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


use = [(1:8) (10:18)]; % Exclude VLF, stdResp, second stdRR, meanResp


%    
%% Regression on all participants for MIST

% mist_stress = zeros(9,5);

mist_stress = horzcat([ 1 3 7 7]',...
               [ 1 6 5 7]',...
               [ 3 5 6 7]',...
               [ 3 5 5 5]',...
               [ 7 8 9 9]',...
               [ 6 8 5 8]',...
               [ 6 8 8 7]',...
               [ 6 8 9 9]');

   
rest=[p1(1,use); ...
    p2(1,use); ...
    p3(1,use); ...
    p4(1,use); ...
    p6(1,use); ...
    p7(1,use); ...
    p8(1,use); ...
    p9(1,use)];
           
X= [p1(2:5,use)./repmat(rest(1,:),4,1); ...
    p2(2:5,use)./repmat(rest(2,:),4,1); ...
    p3(2:5,use)./repmat(rest(3,:),4,1); ...
    p4(2:5,use)./repmat(rest(4,:),4,1); ...
    p6(2:5,use)./repmat(rest(5,:),4,1); ...
    p7(2:5,use)./repmat(rest(6,:),4,1); ...
    p8(2:5,use)./repmat(rest(7,:),4,1); ...
    p9(2:5,use)./repmat(rest(8,:),4,1)];

    
Y = mist_stress; % exclude rest state
Y = mist_stress(2:end,:)./repmat(mist_stress(1,:),3,1); % exclude rest state

xvu = xvars(use);

model = cjw_pls2(X, Y, 10);

%%
vip_all = zeros(size(X,1)/4, size(X,2));
b_coeff_all = zeros(size(X,1)/4, size(X,2));
LV1_p_all = zeros(size(X,1)/4, size(X,2));
R2Y_all = [];
for i=1:4:size(X,1)
    ii = ceil(i/4);
    this_X = [X(i,:); X(i+1,:); X(i+2,:); X(i+3,:)];
    this_Y = mist_stress(:,ii);
    this_model = cjw_pls2(this_X, this_Y, 4);
    this_vip = cjw_vip(this_X, this_Y, this_model);
    this_LV1_p = this_model(1).p';
    LV1_p_all(ii,:) = this_LV1_p;
    R2Y_all(ii) = this_model(1).R2Y;
    b_coeff_all(ii, :) = this_model(1).b_coeff;
    vip_all(ii, :) = this_vip;
    this_vip_lv1 = cjw_vip(this_X, this_Y, this_model);
    vip_lv1_all(ii, :) = this_vip_lv1;
end

%% Boxplot of b_pls coefficients
boxplot(b_coeff_all);
hold on;
axis manual;
plot([0 20], 0*[mean(mean(abs(b_coeff_all))) mean(mean(abs(b_coeff_all)))], 'g');
set(gca, 'XTickLabels', xvu, 'XTickLabelRotation', 90);
ylabel('$b_\textrm{PLS}$');
% bar([abs(mean(b_coeff_all)); mean(vip_all)]', 'group');
% legend('$b_\textrm{PLS}$', 'VIP', 'location', 'eastoutside');
% CJWPlotV4('MIST_Normalised_All_Participant_Model_b_box' );  

%% Boxplot of VIP
boxplot(vip_all);
hold on;
axis manual;
plot([0 20], [1 1], 'g');
ylabel('VIP');
set(gca, 'XTickLabels', xvu, 'XTickLabelRotation', 90);
hold off;
ylim([0 2.2]);
% CJWPlotV4('MIST_Normalised_All_Participant_Model_vip_box' );  

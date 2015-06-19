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

X= [p1(2:5,use); ...
    p2(2:5,use); ...
    p3(2:5,use); ...
    p4(2:5,use); ...
    p6(2:5,use); ...
    p7(2:5,use); ...
    p8(2:5,use); ...
    p9(2:5,use)];
   
    
Y = mist_stress; % exclude rest state

xvu = xvars(use);

vip_all = zeros(size(X,1)/4, size(X,2));
b_coeff_all = zeros(size(X,1)/4, size(X,2));
LV1_p_all = zeros(size(X,1)/4, size(X,2));
R2Y_all = [];
for i=1:4:size(X,1)
    ii = ceil(i/4);
    this_X = [X(i,:); X(i+1,:); X(i+2,:); X(i+3,:)];
    this_Y = mist_stress(:,ii);
    this_model = cjw_pls(this_X, this_Y, 4);
    this_vip = cjw_vip(this_X, this_Y, this_model);
    this_LV1_p = this_model(1).p';
    LV1_p_all(ii,:) = this_LV1_p;
    R2Y_all(ii) = this_model(1).R2Y;
    b_coeff_all(ii, :) = this_model(1).b_coeff;
    b_all(ii) = this_model(1).b;
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
% CJWPlotV4('MIST_All_Participant_Model_b_box' );  

%% Boxplot of VIP
boxplot(vip_all);
hold on;
axis manual;
plot([0 20], [1 1], 'g');
ylabel('VIP');
set(gca, 'XTickLabels', xvu, 'XTickLabelRotation', 90);
hold off;
ylim([0 2.2]);
% CJWPlotV4('MIST_All_Participant_Model_vip_box_DNM' );  

%% Weighted VIPs by R^2... not actually used.
% vip_all_w = zeros(size(vip_all));
% for j=1:size(vip_all,1)
%     vip_all_w(j,:) = vip_all(j,:)*R2Y_all(j);
% end

%% Bar chart of R^2_Y for each of the participant models.
bar(R2Y_all*100);
xlabel('Participant');
ylabel('Response Variance Explained (\%)');
ylim([0 110]);
axis manual;
hold on;
plot([-10 110], [90 90], 'g');
hold off;
set(gca, 'YTickLabel', [arrayfun(@num2str, 0:10:100, 'uniform', false), {''}] );
set(gca, 'XTickLabel', arrayfun(@num2str, [1:4 6:9], 'uniform', false) );
% set(gca, 'XTickLabels', xvu, 'XTickLabelRotation', 90);
% sumTotals = sum(data,2);
% sumTotals = sumTotals(~isnan(sumTotals));
labels = arrayfun(@(x) sprintf('%.2f\\%%', x*100),R2Y_all,'uniform',false);
hText = text(1:size(R2Y_all,2), R2Y_all*100, labels);
set(hText, 'VerticalAlignment','bottom', 'Color','k', 'HorizontalAlignment', 'center');
% CJWPlotV4('MIST_All_Participant_Model_LV1_R2Y_Bar' );  

%% Boxplot of LV1
boxplot(LV1_p_all);
hold on;
plot([0 20], [0 0], 'g');
set(gca, 'XTickLabels', xvu, 'XTickLabelRotation', 90);

% %% Remodel using only participants where R^2_Y of LV1 > 0.9
% 
% nreduced= sum(R2Y_all > 0.9);
% 
% vip_all_r = zeros(nreduced, size(X,2));
% b_coeff_all_r = zeros(nreduced, size(X,2));
% LV1_p_all_r = zeros(nreduced, size(X,2));
% R2Y_all_r = [];
% j=0;
% for i=find(R2Y_all > 0.9)
%     start = (i-1)*4+1;
%     j = j+1;
%     this_X = [X(start,:); X(start+1,:); X(start+2,:); X(start+3,:)];
%     this_Y = mist_stress(:,i);
%     this_model_r = cjw_pls2(this_X, this_Y, 1);
%     this_vip_r = cjw_vip(this_X, this_Y, this_model_r);
%     this_LV1_p_r = this_model_r(1).p';
%     LV1_p_all_r(j,:) = this_LV1_p_r;
%     R2Y_all_r(j) = this_model_r(1).R2Y;
%     b_coeff_all_r(j, :) = this_model_r(1).b_coeff;
%     vip_all_r(j, :) = this_vip_r;
%     this_vip_lv1_r = cjw_vip(this_X, this_Y, this_model_r);
%     vip_lv1_all_r(j, :) = this_vip_lv1_r;
% end
% 

%% Look at the coefficients only for the participants where LV1 explains>90%
figure;
boxplot(b_coeff_all(R2Y_all>0.9,:));
hold on;
axis manual;
plot([0 20], [0 0], 'g');
ylim(0.6*[-1 1]);
ylabel('$b_\textrm{PLS}$');
set(gca, 'XTickLabels', xvu, 'XTickLabelRotation', 90);
hold off;
% ylim([0 2.2]);
% CJWPlotV4('MIST_Some_Participant_Model_LV1_b_box' );  


%% Look at the VIP only for the participants where LV1 explains>90%
figure;
boxplot(vip_all(R2Y_all>0.9,:));
hold on;
axis manual;
plot([0 20], [1 1], 'g');
ylim([0 2.2]);
ylabel('VIP');
set(gca, 'XTickLabels', xvu, 'XTickLabelRotation', 90);
hold off;
% CJWPlotV4('MIST_Some_Participant_Model_LV1_vip_box' );  

%% Reduced model
reduce = find(strcmp(xvars, 'stdHR')|strcmp(xvars, 'LF')|strcmp(xvars, 'LFN')|strcmp(xvars, 'meanIBI'));
Xr=[p1(2:5,use); ...
    p2(2:5,use); ...
    p3(2:5,use); ...
    p4(2:5,use); ...
    p6(2:5,use); ...
    p7(2:5,use); ...
    p8(2:5,use); ...
    p9(2:5,use)];
xvr = xvars(reduce);

vip_all = zeros(size(Xr,1)/4, size(Xr,2));
b_coeff_all = zeros(size(Xr,1)/4, size(Xr,2));
LV1_p_all = zeros(size(Xr,1)/4, size(Xr,2));
R2Y_all = [];
for i=1:4:size(Xr,1)
    ii = ceil(i/4);
    this_X = [Xr(i,:); Xr(i+1,:); Xr(i+2,:); Xr(i+3,:)];
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

%%
plot(zscore(X), 'o');
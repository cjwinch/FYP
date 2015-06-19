addpath('/Users/chris/Documents/FYP/MATLAB/cjw_pls_func');
allxvars = {'meanRR', 'stdRR',    ...
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

xvars = allxvars(use);
%    
%% DA on all participants Stroop

stroop_str = zeros(16,2);
stroop_str(1,:) =  [1 0];
stroop_str(2,:) =  [0 1];
stroop_str(3,:) =  [1 0];
stroop_str(4,:) =  [0 1];
stroop_str(5,:) =  [1 0];
stroop_str(6,:) =  [0 1];
stroop_str(7,:) =  [1 0];
stroop_str(8,:) =  [0 1];
stroop_str(9,:) =  [1 0];
stroop_str(10,:) = [0 1];
stroop_str(11,:) = [1 0];
stroop_str(12,:) = [0 1];
stroop_str(13,:) = [1 0];
stroop_str(14,:) = [0 1];
stroop_str(15,:) = [1 0];
stroop_str(16,:) = [0 1];


X= [p1(7:8,use); ...
    p2(7:8,use); ...
    p3(7:8,use); ...
    p4(7:8,use); ...
    p6(7:8,use); ...
    p7(7:8,use); ...
    p8(7:8,use);
    p9(7:8,use);];
   
%% Generate models for every participant
vip_all = zeros(size(X,1)/2, size(X,2));
b_coeff_all = zeros(size(X,1)/2, size(X,2));
LV1_p_all = zeros(size(X,1)/2, size(X,2));
R2Y_all = [];
for i=1:2:size(X,1)
    ii = ceil(i/2);
    this_X = [X(i,:); X(i+1,:);];
    this_Y = stroop_str(i:i+1,:);
    this_model = cjw_pls(this_X, this_Y, 2);
    this_vip = cjw_vip(this_X, this_Y, this_model);
    this_LV1_p = this_model(1).p';
    LV1_p_all(ii,:) = this_LV1_p;
    R2Y_all(ii) = this_model(1).R2Y;
%     b_coeff_all(ii, :) = this_model(1).b_coeff;
    vip_all(ii, :) = this_vip;
    this_vip_lv1 = cjw_vip(this_X, this_Y, this_model);
    vip_lv1_all(ii, :) = this_vip_lv1;
end
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
   
    
Y = stroop_str; % exclude rest state
 
model = cjw_pls(X, Y, 15);

model2 = cjw_pls(X, Y, 2);

% xvars(vvip>0.8)

%%
vip = cjw_vip(X, Y, model);
bar(vip); hold on;
axis manual;
plot([0 length(vip)+1], [1 1], 'g');
hold off;
set(gca, 'XTick', 1:length(vip));
set(gca, 'XTickLabels', xvars, 'XTickLabelRotation', 90);
ylabel('VIP');
% CJWPlotV4('Stroop_Not_Normalised_VIP');

%% Cross Validation For Only 1 LV using a reduced set of variables
[~,avipi] = sort(vip, 'descend');

for k=1:length(avipi)
    Xrs = X(:,avipi(1:k));

    rleaveout = eye(size(Xrs,1));
    rvip_all = zeros(size(Xrs));
    for i=1:size(Xrs,1)
       l = rleaveout(:,i);
       rthis_x = Xrs(~l,:);
       rthis_y = Y(~l,:);
       rthis_model = cjw_pls2(rthis_x, rthis_y,1);
       rthis_vip = cjw_vip(rthis_x, rthis_y, rthis_model);
       rvip_all(i, :) = rthis_vip;
       %rLV1_b(i) = rthis_model.b;
       %rLV1_t(i, :) = [rthis_model.t];
       %rLV1_u(i, :) = [rthis_model.u];
       %rLV1s(i, :) = [rthis_model.p];
       %rLV1s_out(i, :) = [rthis_model.q];

       rtopredict = [1 Xrs(i,:)];
       rpred = [1 Xrs(i,:)]*rthis_model(1).b_coeff_star;
       [~,rguess]=min(abs(rpred-1),[],2);
       ractual = find(Y(i,:));
       rcorrect(i) = (ractual==rguess);
    end

%     rsqy(k) = 
    acc_r(k) = sum(rcorrect)/length(rcorrect);
end

[hAx, hLine1, hLine2] = plotyy(1:length(avipi), vip(avipi), 1:length(acc_r), acc_r, 'plot', 'plot');

xlim(hAx(1), [1 length(acc_r)]);
xlim(hAx(2), [1 length(acc_r)]);

ylabel(hAx(1), 'Minimum VIP of Included Variables');
ylabel(hAx(2), 'Classification Accuracy');
xlabel('Variables Included');

CJWPlotV4('Stroop_Excluding_Variables');

%% Plot the R2X and R2Y against model order
plot(cumsum([model.R2X]));
hold on;
plot(cumsum([model.R2Y]));
hold off;
xlabel('Number of Latent Variables');
legend('$R^2_X$', '$R^2_Y$');
ylim([0 1]);
% CJWPlotV4('Stroop_Not_Normalised_R2');



%% Try using only the ones with vip>1

xvr = xvars(vip>1.05)
Xr = X(:, vip>1);
modelr = cjw_pls2(Xr, Y, 10);

[PRESSr, RESSr, c_rater] = cjw_plscv(Xr,Y, 'DA')


%% Cross Validation of reduced model
leaveout = eye(size(Xr,1));
vip_all = zeros(size(Xr));
for i=1:size(Xr,1)
   l = leaveout(:,i);
   this_x = Xr(~l,:);
   this_y = Y(~l,:);
   this_model = cjw_pls(this_x, this_y,15);
   this_vip = cjw_vip(this_x, this_y, this_model);
%    this_vip = cjw_vip(this_x, this_y, this_model);
   vip_all(i, :) = this_vip;
%    b_coeff_all(i, :) = this_model(1).b_coeff;
   topredict = [1 Xr(i,:)];
   pred = [1 Xr(i,:)]*this_model(1).b_coeff_star;
   [~,guess]=min(abs(pred-1),[],2);
   actual = find(Y(i,:));
   correct(i) = (actual==guess);
end
sum(correct)/length(correct)



%% Cross Validation of full model
leaveout = eye(size(X,1));
vip_all = zeros(size(X));
for i=1:size(X,1)
   l = leaveout(:,i);
   this_x = X(~l,:);
   this_y = Y(~l,:);
   this_model = cjw_pls2(this_x, this_y,15);
   this_vip = cjw_vip(this_x, this_y, this_model);
%    this_vip = cjw_vip(this_x, this_y, this_model);
   vip_all(i, :) = this_vip;
%    b_coeff_all(i, :) = this_model(1).b_coeff;
   topredict = [1 X(i,:)];
   pred = [1 X(i,:)]*this_model(1).b_coeff_star;
   [~,guess]=min(abs(pred-1),[],2);
   actual = find(Y(i,:));
   correct(i) = (actual==guess);
end
sum(correct)/length(correct)


%% Biplot for X %%%%%%!!!!
% Biplot for X

ObsLabels = {'1C', '1IC', '2C','2IC','3C','3IC','4C','4IC','5C','5IC','6C','6IC','7C','7IC', '8C', '8IC'};

p = [model.p];
t = [model.t];
biplot(p(:,1:2), 'scores', t(:,1:2), ...
       'varlabels', xvars, 'ObsLabels', ObsLabels);
   
   xlabel('Latent Variable 1');
xlabel('Latent Variable 2');
   
% %CJWPlotV4('Stroop_Biplot');

%% Biplot for Reduced

% ObsLabels = {'1C', '1IC', '2C','2IC','3C','3IC','4C','4IC','5C','5IC','6C','6IC','7C','7IC', '8C', '8IC'};

pr = [modelr.p];
tr = [modelr.t];
biplot(pr(:,1:2), 'scores', tr(:,1:2), ...
       'varlabels', xvr);%, 'ObsLabels', ObsLabels);
   
xlabel('Latent Variable 1');
xlabel('Latent Variable 2');
   
% %CJWPlotV4('Stroop_Biplot');

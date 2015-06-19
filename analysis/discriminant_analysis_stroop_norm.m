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
   
allall = [p1(1:5,:); p2(1:5,:); p3(1:5,:); p4(1:5,:); zeros(size(p4(1:5,:))); p6(1:5,:);
    p7(1:5,:); p8(1:5,:); p9(1:5,:)];

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
%
% but this time try normalising the stressed state as a ratio of the
% relaxed state
%

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

Xs(1,:) =  X(1,:)./ X(1,:);
Xs(2,:) =  X(2,:)./ X(1,:);
Xs(3,:) =  X(3,:)./ X(3,:);
Xs(4,:) =  X(4,:)./ X(3,:);
Xs(5,:) =  X(5,:)./ X(5,:);
Xs(6,:) =  X(6,:)./ X(5,:);
Xs(7,:) =  X(7,:)./ X(7,:);
Xs(8,:) =  X(8,:)./ X(7,:);
Xs(9,:) =  X(9,:)./ X(9,:);
Xs(10,:) = X(10,:)./X(9,:);
Xs(11,:) = X(11,:)./X(11,:);
Xs(12,:) = X(12,:)./X(11,:);
Xs(13,:) = X(13,:)./X(13,:);
Xs(14,:) = X(14,:)./X(13,:);
Xs(15,:) = X(15,:)./X(15,:);
Xs(16,:) = X(16,:)./X(15,:);

Y = stroop_str; % exclude rest state
 
model = cjw_pls2(Xs, Y, 15);
 

vvip = cjw_vip(Xs, Y, model);


% VIP=vip(X,Y,[model.t],[model.w],[model.q])
xvars(vvip>1)


%%
vip = cjw_vip(X, Y, model2);
bar(vip); hold on;
axis manual;
plot([0 length(vip)+1], [1 1], 'g');
hold off;
set(gca, 'XTick', 1:length(vip));
set(gca, 'XTickLabels', xvars, 'XTickLabelRotation', 90);
ylabel('VIP');

% CJWPlotV4('Stroop_Normalised_VIP');

%% Plot the R2X and R2Y against model order
plot(cumsum([model.R2X]));
hold on;
plot(cumsum([model.R2Y]));
hold off;
xlabel('Number of Latent Variables');
legend('$R^2_X$', '$R^2_Y$');
ylim([0 1]);
% %CJWPlotV4('Stroop_Normalised_R2');


%%
% Biplot for X
CIC = {'C', 'IC'};
ObsLabels = repmat(CIC,1,8);

p = [model.p];
t = [model.t];
h=biplot(p(:,1:2), 'scores', t(:,1:2), ...
       'varlabels', xvars, 'ObsLabels', ObsLabels);

% o=findobj(h, 'Type', 'Line');
% 
% for i=1:length(o)
%     o(i)
%    if strcmp(o(i).Tag,'obsmarker')
%        xc = o(i).XData(1)+0.05;
%        yc = o(i).YData(1);
%        
%        text(xc,yc,ObsLabels(o(i).UserData));
%    end
% end

xlabel('Latent Variable 1');
xlabel('Latent Variable 2');
   
xlim([-4.25 4.25]);

% %CJWPlotV4('Stroop_Normalised_Biplot');
   
%%
% Biplot for Y


q = [model.q];
u = [model.u];
h=biplot(q(:,1:2), ...
       'varlabels', {'C', 'IC'}, 'ObsLabels', ObsLabels);

% o=findobj(h, 'Type', 'Line');
% 
% for i=1:length(o)
%     o(i)
%    if strcmp(o(i).Tag,'obsmarker')
%        xc = o(i).XData(1)+0.05;
%        yc = o(i).YData(1);
%        
%        text(xc,yc,ObsLabels(o(i).UserData));
%    end
% end

xlabel('Latent Variable 1');
ylabel('Latent Variable 2');
   
% xlim([-4.25 4.25]);

% CJWPlotV4('Stroop_Normalised_Biplot_Y');
   
%% Inner Relation
ir1 = model(1).b;
ir2 = model(2).b;
plot([-1 1], ir1*[-1 1]);
axis tight;
axis manual;
hold on;
plot([-1 1], ir2*[-1 1]);
plot([0 0], 5*[-1 1], 'k');
plot([-1 1], [0 0], 'k');
hold off;
grid on;
legend(sprintf('LV1 ($b=%.2f$)', ir1), ...
       sprintf('LV2 ($b=%.2f$)', ir2), ...
       'Location', 'southeast');
xlabel('Predictor LV')
ylabel('Response LV')

% CJWPlotV4('Stroop_Normalised_Inner_Relation');

%% 3D Biplot
p = [model.p];
t = [model.t];
biplot(p(:,1:3), 'scores', t(:,1:3), ...
       'varlabels', xvars);
   
%% Cross Validation FOR FULL MODEL
leaveout = eye(size(Xs,1));
vip_all = zeros(size(Xs));
for i=1:size(Xs,1)
   l = leaveout(:,i);
   this_x = Xs(~l,:);
   this_y = Y(~l,:);
   this_model = cjw_pls2(this_x, this_y,10);
   this_vip = cjw_vip(this_x, this_y, this_model);
   vip_all(i, :) = this_vip;
  
%    b_coeff_all(i, :) = this_model(1).b_coeff;
   topredict = [1 Xs(i,:)];
   pred = [1 Xs(i,:)]*this_model(1).b_coeff_star;
   [~,guess]=min(abs(pred-1),[],2);
   actual = find(Y(i,:));
   correct(i) = (actual==guess);
end
sum(correct)/length(correct)

%%
[ PRESS, RESS, c_rate ] = cjw_plscv(Xs,Y, 'DA');

[hAx, hLine1, hLine2] = plotyy(1:length(PRESS), PRESS, 1:length(PRESS), c_rate, 'plot', 'plot');

xlim(hAx(1), [1 length(PRESS)]);
xlim(hAx(2), [1 length(PRESS)]);
ylim(hAx(1), [1 PRESS(end-1)]);

ylabel(hAx(1), 'PRESS');
ylabel(hAx(2), 'Classification Accuracy');
xlabel('Model Order');

% %CJWPlotV4('Stroop_Normalised_Crossvalidation_DNM');

%% Try using only the ones with vip>1
xvars(vvip>0.83)
Xr = Xs(:, vvip>0.83);
modelr = cjw_pls2(Xr, Y, 10);
[ PRESSr, RESSr, c_rater ] = cjw_plscv(Xr,Y, 'DA');

%% Cross Validation For Only 1 LV
leaveout = eye(size(Xs,1));
vip_all = zeros(size(Xs));
for i=1:size(Xs,1)
   l = leaveout(:,i);
   this_x = Xs(~l,:);
   this_y = Y(~l,:);
   this_model = cjw_pls2(this_x, this_y,1);
   this_vip = cjw_vip(this_x, this_y, this_model);
   vip_all(i, :) = this_vip;
   LV1_b(i) = this_model.b;
   LV1_t(i, :) = [this_model.t];
   LV1_u(i, :) = [this_model.u];
   LV1s(i, :) = [this_model.p];
   LV1s_out(i, :) = [this_model.q];
%    b_coeff_all(i, :) = this_model(1).b_coeff;
   topredict = [1 Xs(i,:)];
   pred = [1 Xs(i,:)]*this_model(1).b_coeff_star;
   [~,guess]=min(abs(pred-1),[],2);
   actual = find(Y(i,:));
   correct(i) = (actual==guess);
end

%% LV1 ONLY for the input block
hb = bar(mean(-LV1s));
hold on;
set(gca, 'XTick', 1:length(LV1s), 'XTickLabels', xvars, 'XTickLabelRotation', 90);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,mean(-LV1s),std(-LV1s),'r.', 'Marker', 'None')
end
axis manual;
plot([-1 length(LV1s)+1], [0 0], 'k');
hold off;

ylabel('Loading on LV1', 'interpreter', 'latex');
xlabel('Predictor Variable', 'interpreter', 'latex');

% %CJWPlotV4('Stroop_Normalised_LV1_X_Loadings');
box off;
grid on;

set(gcf, 'units', 'inches', 'pos', [0 0 9 3])
set(gcf,'color','w');

% export_fig 'Stroop_Normalised_LV1_X_Loadings.pdf' -q101 


%% LV1 ONLY for the output block
hb = bar(mean(-LV1s_out));
set(gca, 'XTick', 1:length(LV1s_out), 'XTickLabels', {'C', 'IC'}, 'XTickLabelRotation', 90);
axis manual
hold on;
plot([-1 length(LV1s_out)+1], [0 0], 'k');
hold off;

ylabel('Loading on Latent Variable 1');
xlabel('Response Variable');

% %CJWPlotV4('Stroop_Normalised_LV1_Y_Loadings_DNM', 'clean', 0);

%% LV1 Inner Relation
scatter(mean(LV1_t), mean(LV1_u), 'r');
l = [-1 1]*max(abs(mean([LV1_t LV1_u])));
hold on;
plot(5*[-1 1], 5*[-1 1]*mean(LV1_b), 'g');
plot(5*[-1 1], 5*[-1 1]*(mean(LV1_b)+std(LV1_b)), 'Color', [0.6 0.2 0.2]);
plot(5*[-1 1], 5*[-1 1]*(mean(LV1_b)-std(LV1_b)), 'Color', [0.2 0.2 0.6]);
hold off

xlabel('LV1 X ($\vec{t}_1$');
ylabel('LV1 Y ($\vec{u}_1$');

xlim(l);
ylim(l);
text(0.2, 0.4, sprintf('$b=%.2f$',mean(LV1_b)), 'Color', 'g');
legend('Average', 'Average $+$ SD', 'Average $-$ SD');
% %CJWPlotV4('Stroop_Normalised_LV1_Inner_Relation_DNM', 'clean', 0);

%% Cross Validation For Only 1 LV using a reduced set of variables
avip = mean(-LV1s);
[~,avipi] = sort(avip, 'descend');

for k=1:length(avipi)
    Xrs = Xs(:,avipi(1:k));

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

[hAx, hLine1, hLine2] = plotyy(1:length(avipi), avip(avipi), 1:length(acc_r), acc_r, 'plot', 'plot');

xlim(hAx(1), [1 length(acc_r)]);
xlim(hAx(2), [1 length(acc_r)]);

ylabel(hAx(1), 'Minimum VIP of Included Variables');
ylabel(hAx(2), 'Classification Accuracy');
xlabel('Variables Included');

% CJWPlotV4('Stroop_Normalised_Exclude_Vars');


%%
% avip = mean(vip_all);
% [~,avipi] = sort(avip, 'descend');
% 
% for k=1:length(avipi)
%     Xrs = Xs(:,avipi(1:k));
% 
%     rmodel = cjw_pls2(Xrs, Y,1);
%     rvip = cjw_vip(rthis_x, rthis_y, rthis_model);
% 
%     rsqy(k) = rmodel(1).R2Y;
% %     acc_r(k) = sum(rcorrect)/length(rcorrect);
% end
% 
% plot(acc_r);
% hold on;
% plot(avip(avipi));
% hold on;
Yvars = {'Hedonic', 'Goes with Meat', 'Goes with Dessert'}
Y=[...
   14 7 8
   10 7 6
    8 5 5
    2 4 7
    6 2 4
    ];

Xvars = {'Price', 'Sugar', 'Alcohol', 'Acidity'}
X=[...
    7  7 13 7
    4  3 14 7
   10  5 12 5
   16  7 11 3
   13  3 10 3
    ];

wine_n = @(n) sprintf('Wine %d',n);
wine_labels = arrayfun(wine_n, 1:length(Y), 'UniformOutput', false);

model = cjw_pls(X,Y,5);

%% Scree
plot(cumsum([model.R2X]));
hold on;
plot(cumsum([model.R2Y]));
hold off;
ylim([0 1]);
xlabel('Number of Factors');
ylabel('$R^2$');
legend('$R^2_X$', '$R^2_Y$', 'location', 'southeast');
% cjw_press(X,Y);
grid on;
CJWPlotV4('Wine_Scree_Plot');

%% VIP
VIP = cjw_vip(X, Y, model);
bar(VIP);
hold on;
axis manual;
plot([0 5], [1 1], 'g');
hold off;
set(gca, 'XTickLabel', Xvars);
xlabel('Predictor Variable');
ylabel('VIP');
CJWPlotV4('Wine_VIP_Plot');

%% Crossval
[ PRESS, RESS] = cjw_plscv( X, Y );

plot(PRESS);
ylabel('PRESS');
xlabel('Number of Factors');
set(gca, 'XTick', 1:4);
CJWPlotV4('Wine_PRESS_Plot');

%% Biplot for X
figure(1)
p = [model.p];
t = [model.t];
biplot(p(:,1:2), 'scores', t(:,1:2), ...
       'varlabels', Xvars, 'ObsLabels', wine_labels);

xlabel('Latent Variable 1');
ylabel('Latent Variable 2');

% CJWPlotV4('Wine_X_Biplot_DONTMODIFY')

%% Biplot for Y
figure(2)
q = [model.q];
u = [model.u];
biplot(q(:,1:2), 'scores', u(:,1:2), ...
       'varlabels', Yvars, 'ObsLabels', wine_labels);

xlabel('Latent Variable 1');
ylabel('Latent Variable 2');

% CJWPlotV4('Wine_Y_Biplot_DONTMODIFY')

%% Biplot for the inner relationship
q = [model.q];
u = [model.u];

p = [model.p];
t = [model.t];

% tu = [t(:,1) u(:,1)];
% qp = [p(:,1) [q(:,1)]];



figure(1);
scatter(t(:,1), u(:,1), 'r');
l = [-1 1]*max(abs([t(:,1); u(:,1)]));
hold on;
plot(5*[-1 1], 5*[-1 1]*model(1).b, 'g');
xlim(l);
ylim(l);
hold off;
xlabel('$t_1$');
ylabel('$u_1$');
h = legend('Sample data in LV1 space', sprintf('$b=%.2f$',model(1).b), 'Location', 'northwest');
set(h, 'interpreter', 'latex');

% CJWPlotV4('Wine_t_u_1', 'clean', 0)

figure(2);
scatter(t(:,2), u(:,2), 'r');
l = [-1 1]*max(abs([t(:,2); u(:,2)]));
hold on;
plot(5*[-1 1], 5*[-1 1]*model(2).b, 'g');
hold off;
xlim(l);
ylim(l);
xlabel('$t_2$');
ylabel('$u_2$');
h = legend('Sample data in LV2 space', sprintf('$b=%.2f$',model(2).b), 'Location', 'northwest');
set(h, 'interpreter', 'latex');
% CJWPlotV4('Wine_t_u_2', 'clean', 0)

%%

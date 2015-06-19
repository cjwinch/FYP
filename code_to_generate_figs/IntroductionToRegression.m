%% Univariate Regression
% y = bx
clear all

x = 1:20;
y = 2*x + 1.5*randn(size(x));

b = regress(y',x');

scatter(x,y, 'b');
hold on;
plot([0 22], b*[0 22], 'r');
xlim([0 22]);
ylim([0 45]);
hold off;

xlabel('Dependent Variable (x)');
ylabel('Independent Variable (y)');

CJWPlotV4('Univariate_Regression_Example')

%% Multivariate Regression
load hald.mat

% Low col
figure(1);
scatter3(ingredients(:,3), ingredients(:,4), heat)
lim1 = axis;


figure(2);
scatter3(ingredients(:,2), ingredients(:,4), heat)
lim2 = axis;

l = [min(lim1(1), lim2(1)) max(lim1(2), lim2(2)) min(lim1(3), lim2(3))...
     max(lim1(4), lim2(4)) min(lim1(5), lim2(5)) max(lim1(6), lim2(6))];

%% Figure 2
axis(l);
xlabel('3CaO.SiO$_2$ (\%)');
ylabel('2CaO.SiO$_2$ (\%)');
zlabel('Heat (cal/gm)');

b = regress(heat,[ones(size(heat)) ingredients(:,2) ingredients(:,4)]);

x1fit = min(ingredients(:,2)):1:max(ingredients(:,2));
x2fit = min(ingredients(:,4)):1:max(ingredients(:,4));

[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;

hold on;
mesh(X1FIT,X2FIT,YFIT)
colormap winter
hold off;
% CJWPlotV4('Multivariate_Regression_Multicollinear')

%% Figure 1
figure(1);
axis(l);
xlabel('4CaO.Al$_2$O$_3$.Fe$_2$O$_3$ (\%)');
ylabel('2CaO.SiO$_2$ (\%)');
zlabel('Heat (cal/gm)');

b = regress(heat,[ones(size(heat)) ingredients(:,3) ingredients(:,4)]);

x1fit = min(ingredients(:,3)):1:max(ingredients(:,2));
x2fit = min(ingredients(:,3)):1:max(ingredients(:,4));

[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;

hold on;
mesh(X1FIT,X2FIT,YFIT)
colormap winter;
hold off;

view([25.5000, 32]);

% CJWPlotV4('Multivariate_Regression_Non_Multicollinear')

%% Try again with some noise
% Low col

noisy_ingredients = zeros(size(ingredients))
for i=1:size(ingredients,2)
   noisy_ingredients(:,i) = ingredients(:,i) + 0.01*var(ingredients(:,i))*randn(size(ingredients(:,i)));
end

figure(1);
scatter3(noisy_ingredients(:,3), noisy_ingredients(:,4), heat)
lim1 = axis;


figure(2);
scatter3(noisy_ingredients(:,2), noisy_ingredients(:,4), heat)
lim2 = axis;

l = [min(lim1(1), lim2(1)) max(lim1(2), lim2(2)) min(lim1(3), lim2(3))...
     max(lim1(4), lim2(4)) min(lim1(5), lim2(5)) max(lim1(6), lim2(6))];


axis(l);
xlabel('3CaO.SiO$_2$ (\%)');
ylabel('2CaO.SiO$_2$ (\%)');
zlabel('Heat (cal/gm)');

b = regress(heat,[ones(size(heat)) noisy_ingredients(:,2) noisy_ingredients(:,4)]);

x1fit = min(noisy_ingredients(:,2)):1:max(noisy_ingredients(:,2));
x2fit = min(noisy_ingredients(:,4)):1:max(noisy_ingredients(:,4));

[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;

hold on;
h= surf(X1FIT,X2FIT,YFIT);
set(h, 'edgecolor','none');
colormap winter
view([9.5000 32]);

% CJWPlotV4('Multivariate_Regression_Multicollinear_Noise')

figure(1);
axis(l);
xlabel('4CaO.Al$_2$O$_3$.Fe$_2$O$_3$ (\%)');
ylabel('2CaO.SiO$_2$ (\%)');
zlabel('Heat (cal/gm)');

b = regress(heat,[ones(size(heat)) noisy_ingredients(:,3) noisy_ingredients(:,4)]);

x1fit = min(noisy_ingredients(:,3)):1:max(noisy_ingredients(:,2));
x2fit = min(noisy_ingredients(:,3)):1:max(noisy_ingredients(:,4));

[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;

hold on;
h= surf(X1FIT,X2FIT,YFIT);
set(h, 'edgecolor','none');
colormap winter;

view([25.5000, 42]);


%CJWPlotV4('Multivariate_Regression_Non_Multicollinear_Noise')


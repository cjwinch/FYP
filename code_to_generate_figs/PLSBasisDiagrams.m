%% Some data with various covariance matrices

circle = @(r, n) pol2cart(linspace(0,2*pi,n),ones(1,n)*r);

mu = [0,0];
sigma = [1 0;0 1];
rng default  % For reproducibility
r = mvnrnd(mu,sigma,500);

figure;
plot(r(:,1),r(:,2),'x');
hold on;
[X, Y] = circle(1,100);
plot(X,Y);

axis(max(abs(axis)) *  [-1 1 -1 1]);

%%
mu = [0,0];
sigma = [5 0;0 1];
rng default  % For reproducibility
r = mvnrnd(mu,sigma,500);

c2 = [X' Y'] * sigma;

figure;
plot(r(:,1),r(:,2),'x');
hold on;
plot(c2(:,1), c2(:,2));

axis(max(abs(axis)) *  [-1 1 -1 1]);

%%
mu = [0,0];
scale = [1 0; 0 1];
rotate = [cos(pi/4) sin(pi/4); -sin(pi/4) cos(pi/4)];
sigma = scale * rotate;
sigma = [1 0.5;0.5 1];
rng default;  % For reproducibility
r = mvnrnd(mu,sigma,500);

c3 = [X' Y'] * sigma;

figure(3);
plot(r(:,1),r(:,2),'x');
hold on;
plot(c3(:,1), c3(:,2));


axis(max(abs(axis)) *  [-1 1 -1 1]);

%% The process of transforming white data into sigma data
sigma = [1 0.5;0.5 1];
mu = [0,0];
[V,D] = eig(sigma);

% plot(shape(:,1), shape(:,2));
% hold on;
% plot(trans(:,1), trans(:,2));

% r = mvnrnd(mu,sigma,500);

c4 = [X' Y']*sigma;
c3 = [X' Y'];
% trans = c3*V;
trans = c3*D;
trans2 = trans*V;

% plot(r(:,1),r(:,2),'x'); hold on;
axis(max(abs(axis)) *  [-1 1 -1 1]);
hold on
plot(c4(:,1), c4(:,2));  
% plot(c3(:,1), c3(:,2));  
% plot(trans(:,1), trans(:,2));
plot(trans2(:,1), trans2(:,2)); 
hold off;

%% Perfect circular data diagram:
[X,Y] = circle(1,100);
plot(X,Y,'x');
xlabel('x');
ylabel('y');
grid on
box off
axis(3*[-1 1 -1 1]);

CJWPlotV4('Illustration_of_circular_data');

%% Not-circular data diagram:
circle = @(r, n) pol2cart(linspace(0,2*pi,n),ones(1,n)*r);

% Create a circle of data
[X,Y] = circle(1,100);

% Distort the data
sigma = [2 0.5; 0.5 1];
s = [X' Y'] * sigma;

% Find the distortion in the data
[V,D] = eig(sigma);

% [u,s,v]=svd(sigma);


% Use this to convert the circle back.
[v, d] = eig(sigma);
b = [X' Y'] * sqrt(d) *  v ;

after_scale = [X' Y'] * D ;
after_rotate = after_scale * V';

figure(1)
plot(s(:,1), s(:,2), 'xb');
hold on
plot(after_rotate(:,1), after_rotate(:,2),'xr');
%%
figure(1);
plot(s(:,1), s(:,2), 'xb');
xlabel('x');
ylabel('y');
axis(3*[-1 1 -1 1]);
grid on

figure(2);
plot(s(:,1), s(:,2), 'xb');
xlabel('x');
ylabel('y');
axis(3*[-1 1 -1 1]);
grid on
hold on;
scale = 0;
z = zeros(1,2);
sV = V*D*-eye(2);

% plot(b(:,1), b(:,2),'r')

plot([z;sV(1,:)],[z;sV(2,:)],'g')
% quiver(z,z,sV(1,:),sV(2,:),0,'g')
hold off;
CJWPlotV4('Illustration_of_noncircular_data');

figure(3);
projection = [X' Y'] * d;
plot(projection(:,1), projection(:,2))
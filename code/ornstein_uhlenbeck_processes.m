clear all; close all; clc;
subplot(3,2,1)

%%
% Part (a)
T = 2;
N=100*T; % Number of changes in track
h=sqrt(T/N);
x = zeros(N,1);
o = zeros(N,1);
% Iterate over N timesteps
for i=1:N
    r = randn();
    x(i+1,1)=x(i,1)+h*r;
    o(i+1,1)=o(i,1)+h*r-o(i,1)*.05;
end
% Plot
hold on; 
plot(x);
plot(o);
plot([0,200],[0,0]);
title("Brownian and Ornstein-Uhlenbeck");
legend("Brownian", "Ornstein-Uhlenbeck");
set(gca,'XTick',[],'YTick',[])
hold off;

%%
% Part (b)
T= 10;
N=100*T;
h=sqrt(T/N); 
x = zeros(N,1);
y = zeros(N,1);
% Iteration to store positions of particles
for i=1:N
    x(i+1)=x(i)+h*randn();
    y(i+1)=y(i)+h*randn();
end
% Plot
subplot(3,2,2)
plot(x,y)
title("X, Y Brownian");
set(gca,'XTick',[],'YTick',[])

%%
% Part (c)
n = 10000;
% State variables: x, y, theta
state = zeros(n,3);
for i=2:n
    state(i,3) = state(i-1,3) + randn() * 3; % theta
    state(i,1) = state(i-1,1) + cosd(state(i,3)); % x 
    state(i,2) = state(i-1,2) + sind(state(i,3));% y 
end
subplot(3,2,3)
plot(state(:,1), state(:,2))
title("Theta Brownian");
set(gca,'XTick',[],'YTick',[])

%%
% Part (d)
n = 3000;
% State variables: x, y, dx, dy
state = zeros(n,4); 
offset = 0.01;
for i=2:n
    state(i,3) = state(i-1,3) + randn() - state(i-1,3)*offset; % dx
    state(i,4) = state(i-1,4) + randn() - state(i-1,4)*offset; % dy
    state(i,1) = state(i-1,1) + state(i,3); % x 
    state(i,2) = state(i-1,2) + state(i,4);% y 
end
subplot(3,2,4)
plot(state(:,1), state(:,2))
title("V_x, Y_x Ornstein-Uhlenbeck");
set(gca,'XTick',[],'YTick',[]);
 
%%
% Part (e)
% change_in_deg = len*k;
n = 3000;
len = 0.6;
% x, y, theta, curvature(k)
state = zeros(n,4); 
for i=2:n
    state(i,4) = state(i-1,4) + randn()/3; % curvature
    state(i,3) = state(i-1,3) + len * state(i,4); % theta
    state(i,1) = state(i-1,1) + cosd(state(i,3)); % x 
    state(i,2) = state(i-1,2) + sind(state(i,3));% y 
end
subplot(3,2,5)
plot(state(:,1), state(:,2))
title("Curvature Brownian");
set(gca,'XTick',[],'YTick',[])
 
%%
% Part (f)
n = 4000;
len = 0.6;
% x, y, theta, curvature(k)
state = zeros(n,4); 
for i=2:n
    state(i,4) = state(i-1,4) + (-state(i-1,4)*.05) + randn()/3; % Curvature(k)
    state(i,3) = state(i-1,3) + len * state(i,4); % theta
    state(i,1) = state(i-1,1) + cosd(state(i,3)); % x 
    state(i,2) = state(i-1,2) + sind(state(i,3));% y 
end
subplot(3,2,6)
hold off;
plot(state(:,1), state(:,2))
title("Curvature Ornstein-Uhlenbeck");
set(gca,'XTick',[],'YTick',[])

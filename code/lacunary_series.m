clear all; close all; clc

% Hyperparameters
n = 100; % Choose a reasonable n or else everything overflows
Xes = 100;
x_range = 0:5/Xes:5; % X is (0,5)
a = 0.1:0.2:0.9; % 0 < a < 1
b = 2.1:0.2:2.9; % 1 < b
g = zeros(Xes, 5);

for k = 1:numel(x_range) % Iterate over xes while summing to plot each series
    x = x_range(k);
    sum = zeros(1,5);
    for j = 1:n
        sum = (a.^j) .* sin(x.*b.^j) + sum;   
    end
    g(k,:) = sum;
end
plot(x_range, g)
legend("a=0.1 b=2.1","a=0.3 b=2.3","a=0.5 b=2.5","a=0.7 b=2.7","a=0.9 b=2.9"); 
xlabel("x")
ylabel("y")
title("Lacunary Series")
%The Derivative of g simply adds an additional b^n term: 
% a^n * cos(x*b^n) * b^n
 
gderiv = zeros(Xes,5);
 
for k = 1:numel(x_range) % Iterate over xes while summing to plot each series
    x = x_range(k);
    sum = zeros(1,5);    
    for j = 1:n
        sum = a.^j .* cos(x.*b.^j) .* b.^j + sum;
    end
    gderiv(k,:) = sum;
end
figure
% This plot is dominated by the last b because it blows up faster than the
% others
plot(gderiv)
legend("a=0.1 b=2.1","a=0.3 b=2.3","a=0.5 b=2.5","a=0.7 b=2.7","a=0.9 b=2.9");

figure 
for i = 1:5
	subplot(3,2,i)
	plot(g(:,i));
	hold on;
	plot(gderiv(:,i));
    title(sprintf("a=%1.1f b=%1.1f", a(i), b(i)))
    legend("g(x)", "g'(x)")
end
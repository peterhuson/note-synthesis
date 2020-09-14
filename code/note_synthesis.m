close all; clc; 
data = audioread("fivesec.au");

% Hyperparameters
N = 40000;
K = 35; % 30 > 22 meaning leave space for choosing optimal K
n_interval = 40; % interval = 10 runs for a few minutes
n_range = 1:n_interval:40000;
max_p = 40;
p_range = 5:max_p;

C = zeros(40000,36); 
for i = 1:39960
    for j = 1:36
        p = j + 4;
        C(i,j) = data(i)*data(i+p);
    end
end

for i = 39961:40000 % Make sure the last 40 rows of C are populated 
    for j = 1:36
        p = j+4;
        C(i,j) = data(i)*data(i-p);
    end
end
% C is just a cumulative sum of all the individual products
C = cumsum(C); 

amplitude = sum(data.^2);
%%
% a_range = 0.0006:-0.000025:0.00005;

% Found optimal a to be:
a_range = 0.00015;
max_K = zeros(numel(a_range), 1);
for a_index = 1:numel(a_range)
    a = a_range(a_index);

    likelihood = zeros(numel(n_range),K); % Likelihood for (n,k)
    p_argmax = zeros(numel(n_range),K); % Best p for this (n,k)
    n_argmax = ones(numel(n_range),K); % Best n for this (n,k)

    for j = 2:numel(n_range) % Initialize base case with constant nk-1 
        n = n_range(j);
        [likelihood(j,1),p_argmax(j,1)] = max(C(n,:) - C(1,:)); 
    end

    for k = 2:K % Iterate over K Columns
        for j = k:numel(n_range) % Iterate from k to n
            n = n_range(j); 
            % log likelihood sums the cnk - cnk-1
            maximizer = likelihood(1:j-1,k-1) - a*k + C(n,:) - C(1:n_interval:n-1,:); % C(nk,pk) - C(nk-1,pk)
            % This is fast because maximizer is a 2D array of (nxp)
            likelihood(j,k) = max(maximizer(:));
            % I wish matlab would give you the indeces rather than making
            % you find them
            [nk,p] = find(maximizer==likelihood(j,k)); 
            % If there is more than one match, choose the first 
            %TODO: Make this maximization faster 
            p_argmax(j,k) = p(1); 
            n_argmax(j,k) = nk(1);
        end
    end
    [m, max_K(a_index)] = max(likelihood(N/n_interval,:));
end

%%
% Plot the best a
% figure
% plot(a_range, max_K)
% title("max_K vs a")
% xlabel("a")
% ylabel("K")

%%
K = 22;

% Backtracking: Start at 40000, select the optimal n and p
frequency = zeros(K,1);
starting_n_indeces = zeros(K,1);
n_index_counter = N/n_interval;
for k = K:-1:1
    frequency(k) = p_argmax(n_index_counter,k) + 4; 
    starting_n_indeces(k) = n_argmax(n_index_counter, k);
    n_index_counter = starting_n_indeces(k);
end

% Draw Lines on spectrogram
figure
specgram(data,[],8000);
hold on;
start_seconds = starting_n_indeces.*(n_interval/8000);
frequencies = 8000./frequency;
for k = 1:K
    x = start_seconds(k);
    line([x, x],[0,4000]); % draw vertical lines
    y = frequencies(k);
    if k < K % draw horizontal lines
        line([start_seconds(k) start_seconds(k+1)], [y, y], 'LineWidth', 4); 
    else
        line([start_seconds(k) 5], [y, y], 'LineWidth', 4);
    end
end
zl = zlim;
title(sprintf("Spectrogram K=%d",K));
axis([xlim ylim zl(1) max(0, zl(2))]);
view(0,90);

% Print Optimal Sequence
disp("Boundary(ms)  Frequency(Hz)")
fprintf('%8.0f %8.0f\n', [starting_n_indeces*n_interval/8,frequencies]')
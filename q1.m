% Load all data
clear
load 'ps7_data.mat'

NUM_DATA = size(Spikes, 2);
DIMENSION = size(Spikes, 1);

mean_spike = mean(Spikes, 2);
cov_spike = cov(Spikes');

% Plotting Raw spikes
for i=1:NUM_DATA
    plot(Spikes(:,i))    
    hold on
end
xlabel('Time');
ylabel('Voltage');
title('Raw Spike Snippets');


figure;
% Diagonalizing Covariance Matrix
[U, Lambda] = eig(cov_spike);

% Plotting Three Largest Eigen Values
plot(U(:,end), 'r');
hold on
plot(U(:,end-1), 'g');
hold on
plot(U(:,end-2), 'b');
xlabel('Time');
ylabel('Voltage');
legend('PC1', 'PC2', 'PC3');
title('First Three Principal Components');
figure;

% Calculating and plotting Square Rooted Eigen Values
eigen_values = diag(Lambda);
eigen_values_increasing = fliplr(eigen_values')';
plot(fliplr(sqrt(eigen_values_increasing)));
hold on
plot(fliplr(sqrt(eigen_values_increasing)), '*');
xlabel('Component Number');
ylabel('sqrt - eigen value');
title('Square Rooted Eigen Value Spectrum');
figure


% Creating Scatter Plot
PC1 = -U(:,end);
PC2 = U(:,end-1);
mean_repeated = repmat(mean_spike, 1, NUM_DATA);
Z1 = PC1' * (Spikes - mean_repeated);
Z2 = PC2' * (Spikes - mean_repeated);
plot(Z1, Z2, '.k');
xlabel('1st Component Score');
ylabel('2nd Component Score');
title('Scatter plot of Scores From First 2 components');

% Saving Data for Question 2
save 'q1_data.mat';
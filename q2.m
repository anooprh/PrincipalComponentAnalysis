clear;
load 'q1_data.mat'

NUM_FOLDS = 4;

Z = [Z1;Z2];  % 2 Dimensional Representation Of Spikes

REDUCED_DIMENSION = size(Z,1); 
RANGE_OF_GAUSSIANS = 6;
Z_folded = mat2cell(Z, [REDUCED_DIMENSION], repmat(NUM_DATA/NUM_FOLDS, 1, NUM_FOLDS));



likelihood = zeros(1, RANGE_OF_GAUSSIANS);

for k=1:RANGE_OF_GAUSSIANS
    K = k; % Make it 1:8 later
    params.mu = InitParams.mu(:,1:K);
    params.sigma = repmat(InitParams.Sigma, [1,1,K]);
    params.pi = repmat(1/K,1,K);
    local_cross_validation_sum = 0;
    for i=1:NUM_FOLDS
        [mu, sigma, ppi] = func_GMM(params, Z_folded{i});
        for j=1:NUM_DATA
%             Figure out a way to compute N(xn | mu1,mu2, sigma1,sigma2)
%             local_cross_validation_sum = local_cross_validation_sum + ...
%                 logmvnpdf(Z(:,j), mu, sigma, ppi);
        end
    end
    likelihood(k) = local_cross_validation_sum;
end
plot(likelihood);

% Q2B. For plotting the ellipse .. Still have to figure out why it does not
% work for 8 Gaussians
figure
for k=1:RANGE_OF_GAUSSIANS
    K = k;
    subplot(4,2,k);
    plot(Z1, Z2, '.');
    hold on
    params.mu = InitParams.mu(:,1:K);
    params.sigma = repmat(InitParams.Sigma, [1,1,K]);
    params.pi = repmat(1/K,1,K);
    [mu, sigma, ppi] = func_GMM(params, Z_folded{1});
    for i=1:K
        this_point = mu(:,i);
        plot(this_point(1), this_point(2),'k.');
        hold on
        func_plotEllipse(this_point, sigma(:,:,i));
        hold on
    end
    xlabel('X(1)');
    ylabel('X(2)');
    title_str = sprintf('Clustering for %d Gaussians', K);
    title(title_str);
end

% Problem 2c
figure;
K = 3;
recovered_spikes = zeros(DIMENSION, K);
params.mu = InitParams.mu(:,1:K);
params.sigma = repmat(InitParams.Sigma, [1,1,K]);
params.pi = repmat(1/K,1,K);
[mu, sigma, ppi] = func_GMM(params, Z_folded{1});

% Xn = U_m * Z_n + mean_spike
U_m = fliplr(U(:,end-1:end));
for i =1:K
    recovered_spikes(:, i) = U_m * mu(:,i) + mean_spike; 
end
plot(recovered_spikes(:, 1), 'r'); hold on
plot(recovered_spikes(:, 2), 'g'); hold on
plot(recovered_spikes(:, 3), 'b'); hold on
xlabel('Time');
ylabel('Voltage');
legend('Cluster Center 1', 'Cluster Center 2', 'Cluster Center 3')
title('High D vectors from 2-D space');

function [ bt ] = boxed_term( xn, params, covvar)
% Computes log(sum over all j, N(xn|mu_j, Sigma_j) * pi_j)
    NUM_CLUSTERS = size(params.mu, 2);

    a = zeros(1, NUM_CLUSTERS);
    max_aj = 0;
    for j=1:NUM_CLUSTERS
        a(j) = logmvnpdf(xn, params.mu(:,j), covvar(:,:,j), params.pi(j));
        if a(j) > max_aj
            max_aj = a(j);
        end
    end
    sum = 0;
    for j=1:NUM_CLUSTERS
        sum = sum + exp(a(j) - max_aj); 
    end
    bt = max_aj + log(sum);
end


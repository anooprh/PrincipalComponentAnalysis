function [logp] = logmvnpdf(x,mu,Sigma, pi_k)
% Computes log(N(xn|mu, Sigma) * pi_k) 

    [D,~] = size(x);
    const = -0.5 * D * log(2*3.1416);

    term1 = -0.5 * ((x - mu)' * (inv(Sigma) * (x - mu)));
    term2 = - 0.5 * logdet(Sigma);    
    logp = (const + term1 + term2) +  log(pi_k);
end

function y = logdet(A)
    y = log(det(A));
end
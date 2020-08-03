function cF = combine_K_mu(K,mu)
m = size(K,3);
n = size(K,1);
cF = zeros(n);
for p =1:m
    cF = cF + K(:,:,p)*mu(p);
end
function [H_normalized, mu, G, obj] = optimal_neighborhood_local_kernel_clustering(K,class_num,lambda,tau,rho)

sample_num = size(K,1);
kernel_num = size(K,3);
mu = ones(kernel_num,1)/kernel_num;
K_mu = combine_K_mu(K,mu);
G = K_mu;
% Calculate Neighborhood of each sample
[NS,NS_num] = genarate_neighborhood(K_mu,tau);%% tau*num

Ai_sum = zeros(sample_num);
for i = 1:sample_num
    Ai_sum(NS(1:NS_num(i),i),NS(1:NS_num(i),i)) = Ai_sum(NS(1:NS_num(i),i),NS(1:NS_num(i),i)) + 1;
end
Ai_sum = Ai_sum./sample_num;

Mi_sum = zeros(kernel_num);
M = zeros(kernel_num);
for p = 1:kernel_num
    for q = 1:kernel_num
        Mi_sum(p,q) = trace(Ai_sum.*K(:,:,p)*K(:,:,q));
        M(p,q) = trace(K(:,:,p)'*K(:,:,q));
    end
end


t = 0;
flag = 1;
while flag
    %% update H with G
    %fprintf('update H with G in t: %d\n', t);
    %tic;
    H = my_kernel_kmeans(Ai_sum.*G,class_num);
    %toc;
    %% update mu with H and G
    %fprintf('update mu with H and G in t: %d\n', t);
    %tic;
    [mu] = update_mu(K,G,M,Mi_sum,rho,lambda);
    %toc;
    K_mu = combine_K_mu(K, mu);
    %% update G with K and H
    %fprintf('update G with K and H in t: %d\n', t);
    %tic;
    G = update_G(K_mu,Ai_sum,rho, H);
    %toc;
    t = t+1;
    %fprintf('cal objective t: %d\n', t);
    %tic;
    [obj(t)] = cal_obj(G,H,K_mu,Ai_sum,Mi_sum,mu,rho,lambda);
    %toc;
    if t>=2 && (abs((obj(t-1)-obj(t))/(obj(t-1)))<1e-4 || t>100)
        flag =0;
    end
end
H_normalized = H./ repmat(sqrt(sum(H.^2, 2)), 1,class_num);

end

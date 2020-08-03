% src_path = '/home/ftp2/jiyuan/Kernel_Kmeans/dataset/';
src_path = 'D:\Work\lab_xwl\code\dataset\';
data_name = 'proteinFold';
% data_name = 'flower17';
% data_name = 'bbcsport2view';
% data_name = 'UCI_DIGIT';
% data_name = 'cornell';
load([src_path, data_name, '_Kmatrix'], 'KH', 'Y');

class_num = length(unique(Y));
kernel_num = size(KH, 3);
sample_num = size(KH, 1);

KH = kcenter(KH);
KH = knorm(KH);

lambda = 2;
tau_set = [-0.5:0.1:0.5];
rho_set = 2.^[-1];

acc = zeros(length(tau_set), length(rho_set));
nmi = zeros(length(tau_set), length(rho_set));
pur = zeros(length(tau_set), length(rho_set));

n = 0;

for it = 1:length(tau_set)
    mu = ones(kernel_num,1)/kernel_num;
    K_mu = combine_K_mu(KH,mu);
    % Calculate Neighborhood of each sample
    [NS,NS_num] = genarate_neighborhood(K_mu,tau_set(it));%% tau*num

    Ai_sum = zeros(sample_num);
    for i = 1:sample_num
        Ai_sum(NS(1:NS_num(i),i),NS(1:NS_num(i),i)) = Ai_sum(NS(1:NS_num(i),i),NS(1:NS_num(i),i)) + 1;
    end
    Ai_sum = Ai_sum./sample_num;

    Mi_sum = zeros(kernel_num);
    M = zeros(kernel_num);
    for p = 1:kernel_num
        for q = 1:kernel_num
            Mi_sum(p,q) = trace(Ai_sum.*KH(:,:,p)*KH(:,:,q));
            M(p,q) = trace(KH(:,:,p)'*KH(:,:,q));
        end
    end
    for ir=1:length(rho_set)
        tic;
        [H_normalized, mu, G, obj] = optimal_neighborhood_local_kernel_clustering(KH, class_num, lambda, tau_set(it), rho_set(ir),Y, Ai_sum, Mi_sum, M);
        eval = my_eval(H_normalized, Y, class_num);
        toc;
        acc(it, ir) = eval(1);
        nmi(it, ir) = eval(2);
        pur(it, ir) = eval(3);
        n=n+1;
        fprintf('\nlambda: %d, tau: %d, rho: %d',lambda, tau_set(it), rho_set(ir));
        fprintf('\nacc: %f, nmi: %f, pur: %f\n', eval(1), eval(2), eval(3)); 
    end
end

res = [max(max(max(acc))), max(max(max(nmi))), max(max(max(pur)))];
% res_path = '/home/ftp2/jiyuan/Kernel_Kmeans/ON-ALKA/res/';
res_path = 'D:\Work\work2015\jiyuan\code\res\';
save([res_path,data_name,'_cluster_res_selparam1.mat'],'acc','nmi','pur');

opttau_indx = find(acc==max(acc));
tau_set = [tau_set(opttau_indx)];
rho_set = 2.^[-15:1:15];

acc = zeros(length(tau_set), length(rho_set));
nmi = zeros(length(tau_set), length(rho_set));
pur = zeros(length(tau_set), length(rho_set));

for it = 1:length(tau_set)
    
    mu = ones(kernel_num,1)/kernel_num;
    K_mu = combine_K_mu(KH,mu);
    % Calculate Neighborhood of each sample
    [NS,NS_num] = genarate_neighborhood(K_mu,tau_set(it));%% tau*num

    Ai_sum = zeros(sample_num);
    for i = 1:sample_num
        Ai_sum(NS(1:NS_num(i),i),NS(1:NS_num(i),i)) = Ai_sum(NS(1:NS_num(i),i),NS(1:NS_num(i),i)) + 1;
    end
    Ai_sum = Ai_sum./sample_num;

    Mi_sum = zeros(kernel_num);
    M = zeros(kernel_num);
    for p = 1:kernel_num
        for q = 1:kernel_num
            Mi_sum(p,q) = trace(Ai_sum.*KH(:,:,p)*KH(:,:,q));
            M(p,q) = trace(KH(:,:,p)'*KH(:,:,q));
        end
    end
    for ir=1:length(rho_set)
        tic;
        [H_normalized, mu, G, obj] = optimal_neighborhood_local_kernel_clustering(KH, class_num, lambda, tau_set(it), rho_set(ir),Y, Ai_sum, Mi_sum, M);
        eval = my_eval(H_normalized, Y, class_num);
        toc;
        acc(it, ir) = eval(1);
        nmi(it, ir) = eval(2);
        pur(it, ir) = eval(3);
        fprintf('\nlambda: %d, tau: %d, rho: %d',lambda, tau_set(it), rho_set(ir));
        fprintf('\nacc: %f, nmi: %f, pur: %f\n', eval(1), eval(2), eval(3)); 
    end
end

res2 = [max(max(max(acc))), max(max(max(nmi))), max(max(max(pur)))];
% res_path = '/home/ftp2/jiyuan/Kernel_Kmeans/ON-ALKA/res/';
res_path = 'D:\Work\work2015\jiyuan\code\res\';
save([res_path,data_name,'_cluster_res_selparam2.mat'],'acc','nmi','pur');

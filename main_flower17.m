src_path = '/home/ftp2/jiyuan/Kernel_Kmeans/dataset/';
% src_path = 'D:\Work\code\dataset\';
data_name = 'flower17';
% data_name = 'flower17';
% data_name = 'bbcsport2view';
% data_name = 'UCI_DIGIT';
% data_name = 'cornell';
load([src_path, data_name, '_Kmatrix'], 'KH', 'Y');

class_num = length(unique(Y));
base_kernel_num = size(KH, 3);
sample_num = size(KH, 1);

K = kcenter(KH);
K = knorm(K);

lambda = 2;
tau_set = [-0.5:0.1:0.5];
rho_set = 2.^[-15:1:15];

acc = zeros(length(tau_set), length(rho_set));
nmi = zeros(length(tau_set), length(rho_set));
pur = zeros(length(tau_set), length(rho_set));

n = 0;

for it = 1:length(tau_set)
    for ir=1:length(rho_set)
        tic;
        [H_normalized, mu, G, obj] = optimal_neighborhood_local_kernel_clustering(K, class_num, lambda, tau_set(it), rho_set(ir),Y);
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
res_path = '/home/ftp2/jiyuan/Kernel_Kmeans/ON-ALKA/res/';
save([res_path,data_name,'_cluster_res.mat'],'acc','nmi','pur');
function [ mu,fval,exitflag ] = update_mu( K,G,M,Mi_sum,rho,lambda )

sample_num = size(K,1);
kernel_num = size(K,3);    

a =zeros(kernel_num,1);
for p =1:kernel_num
    a(p) = -rho*trace(G'*K(:,:,p));
end

coff = (lambda*Mi_sum + rho*M);
H = (coff+coff')/2;
A = [];
b = [];
Aeq = ones(1,kernel_num);
beq = 1;
lb = zeros(kernel_num,1);
ub = ones(kernel_num,1);
[mu,fval,exitflag] = quadprog(H,a,A,b,Aeq,beq,lb,ub);
mu = mu./sum(mu);

end


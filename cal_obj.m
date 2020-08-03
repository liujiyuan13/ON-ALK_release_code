function [obj] = cal_obj(G,H,K_mu,Ai_sum,Mi_sum,mu,rho,lambda )

sample_num = size(K_mu,1);

term1 = trace(Ai_sum.*G*(eye(sample_num,sample_num)-H*H'));

term2 = mu'*Mi_sum*mu;

term3 = trace((G-K_mu)*(G-K_mu));

obj = term1 + (lambda/2)*term2 + (rho/2)*term3;

end


function G = update_G(K_mu,Ai_sum,rho, H)

sample_num = size(K_mu, 1);

B = K_mu - 1/rho*Ai_sum.*(eye(sample_num,sample_num) - H*H');
B = (B+B')/2;
[V,D] = eig(B);
diagD = diag(D);
diagD(diagD<eps)=0;
G = V*diag(diagD)*V';
G = (G+G')/2;

end
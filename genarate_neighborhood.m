function [indx_0, indx_num]  = genarate_neighborhood(KC,tau)

num = size(KC,1);
KC0 = KC - 10^8*eye(num);
KC = KC0 > tau;
[val,indx] = sort(KC,'descend');
indx_num = sum(val);
indx_0 = indx(1:max(indx_num),:);

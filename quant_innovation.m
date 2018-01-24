
function [nu, S] = quant_innovation(Ez, Ppred, z_a_b, H, R,Cov_zA,QKF_OptK_flag)
nu = Ez - H*z_a_b; %% innovation
if QKF_OptK_flag == 1
S = R + H*Ppred*H' + Cov_zA; %% innovation covariance
else 
S = R + H*Ppred*H'; %% innovation covariance    
end

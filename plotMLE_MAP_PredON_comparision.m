clear all;
close all;
TestedCFG = 'MLE_MAP_QKF_Comp_test';%fitted Gaussian Estimator
TestedMainTitle = 'MLE, MAP and QKF Comparision ';
TestedMainSubTitle = ' - Fixed Observer (Predict-ON, Estimation-ON, CP-Back, CrossingUpdate-On) ';
set(0,'defaultTextInterpreter','latex'); %trying to set the default
fig_ct = 1;
testCase =10;
figure(10);
hold on;grid on;
load mlePred.mat;
N = length(P);
plot(1:N,squeeze(sqrt(P(1,1,1:N))),':b','LineWidth',2,'MarkerSize',12);%Dotted line
load mapPred.mat;
plot(1:N,squeeze(sqrt(P(1,1,1:N))),'--r','LineWidth',2,'MarkerSize',12);%%Dashed line
load qklPred.mat;
plot(1:N,squeeze(sqrt(P(1,1,1:N))),'-.g','LineWidth',2,'MarkerSize',12);%%Dash-dotted line	

xlabel('time-step [s]','Interpreter','Latex');ylabel('$$\sigma [m]$$','Interpreter','Latex');
legend('$$\hat{\sigma}_{mle}$$','$$\hat{\sigma}_{map}$$','$$\hat{\sigma}_{qkl}$$');  set(legend,'Interpreter','latex');
fig_ct = fig_ct+1;

tit1 = title([TestedMainTitle,TestedMainSubTitle]);
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);

figure(11);
hold on;grid on;
load mlePred.mat;
plot(1:N,xhat(1,1:N),':b','LineWidth',2,'MarkerSize',10);%Dotted line
load mapPred.mat;
plot(1:N,xhat(1,1:N),'--r','LineWidth',2,'MarkerSize',12);%%Dashed line
load qklPred.mat;
plot(1:N,xhat(1,1:N),'-.g','LineWidth',2,'MarkerSize',10);%%Dash-dotted line	

xlabel('time-step [s]','Interpreter','Latex');ylabel('CP Displacment [m]','Interpreter','Latex');
legend('$$\hat{x}_{mle}$$','$$\hat{x}_{map}$$','$$\hat{x}_{qkf}$$');  set(legend,'Interpreter','latex');
fig_ct = fig_ct+1;
tit1 = title([TestedMainTitle,TestedMainSubTitle]);
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);

figure(12);
hold on;grid on;
load mlePred.mat;
plot(1:N,xhat(3,1:N),':b','LineWidth',2,'MarkerSize',10);%Dotted line
load mapPred.mat;
plot(1:N,xhat(3,1:N),'--r','LineWidth',2,'MarkerSize',12);%%Dashed line
load qklPred.mat;
plot(1:N,xhat(3,1:N),'-.g','LineWidth',2,'MarkerSize',10);%%Dash-dotted line	

xlabel('time-step [s]','Interpreter','Latex');ylabel('CP Correlated Spread Rate [m/s]','Interpreter','Latex');
legend('$$\hat{x}_{mle}$$','$$\hat{x}_{map}$$','$$\hat{x}_{qkf}$$');  set(legend,'Interpreter','latex');
fig_ct = fig_ct+1;
tit1 = title([TestedMainTitle,TestedMainSubTitle]);
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);

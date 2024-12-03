
    clear; close all; clc;
%     load('best_ite_5_10.mat')
%     load('best_rmse_5_10.mat')
%     load('bper_ite_5_10.mat')
%     load('bper_rmse_5_10.mat')
    load('2_ite_5_4.mat')
    load('2_rmse_5_4.mat')
    load('ite_pert_5_4.mat')
    load('rmse_pert_5_4.mat')
    figure(2)
%     plot(itebest(:,1),rmsebest(:,1),'b-o','MarkerSize',7,'LineWidth',1)
%     hold on
% 
%     plot(itebper(:,1),rmsebper(:,1),'r--*','MarkerSize',7,'LineWidth',1.5)
    plot(ite2(:,1),rmse2(:,1),'b-o','MarkerSize',7,'LineWidth',1)
    hold on

    plot(ite_pert(:,1),rmse_pert(:,1),'r--*','MarkerSize',7,'LineWidth',1.5)
%     hold on 
% 
%     plot(itebper(:,1),rmsebper(:,1),'g:square','MarkerSize',7,'LineWidth',1.5)
    hold off
    grid on

    title('Tunning T1-PD controller using FA')
    xlabel('Number of Iterations'); 
    ylabel('RMSE for circle tracking');
    legend('Best solutions','Perturbed solutions')
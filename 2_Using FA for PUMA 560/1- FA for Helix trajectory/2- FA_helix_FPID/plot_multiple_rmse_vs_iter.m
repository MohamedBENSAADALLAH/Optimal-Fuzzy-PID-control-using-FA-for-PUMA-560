    clear; close all; clc;
%     load('best_ite_25_8.mat')
%     load('best_rmse_25_8.mat')
%     load('bper_ite_25_8.mat')
%     load('bper_rmse_25_8.mat')
%     load('ite_25_8.mat')
%     load('rmse_25_8.mat')
    load('ite_fin2_25_8.mat')
    load('rmse_fin2_25_8.mat')
    load('ite_pert_25_8.mat')
    load('rmse_pert_25_8.mat')
    
    figure(2)
%     plot(itebest(:,1),rmsebest(:,1),'b-o','MarkerSize',7,'LineWidth',1)
%     hold on
% 
%     plot(itebper(:,1),rmsebper(:,1),'r--*','MarkerSize',7,'LineWidth',1.5)
    plot(ite_fin2(:,1),rmse_fin2(:,1),'b-o','MarkerSize',7,'LineWidth',1)
    hold on

    plot(ite_pert(:,1),rmse_pert(:,1),'r--*','MarkerSize',7,'LineWidth',1.5)
%     hold on 
% 
%     plot(itebper(:,1),rmsebper(:,1),'g:square','MarkerSize',7,'LineWidth',1.5)
    hold off
    grid on

    title('Tunning T1-PID controller using FA')
    xlabel('Number of Iterations'); 
    ylabel('RMSE for helix tracking');
    legend('Best solutions','Perturbed solutions')
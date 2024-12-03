figure

% text = 'RMSE='+str_RMSE+' '+'ISE='+str_ISE
text_metric = Metric(error)
subplot(2,2,1); plot(error); title ('total error');
 legend(text_metric);
% ax2 = gca;
% legend({'y = sin(x)','y = cos(x)'},'Location','southwest')
% hold off; 
% legend(str_RMSE,str_ISE);
hold on;
text_metric='';
text_metric = Metric(error1)
subplot(2,2,2); plot(error1); title ('theta1 error'), hold on; 
 legend(text_metric);
text_metric='';
text_metric = Metric(error2)
subplot(2,2,3); plot(error2); title ('theta2 error'), hold on; 
legend(text_metric);
text_metric='';
text_metric = Metric(error3)
subplot(2,2,4); plot(error3); title ('theta3 error'), hold on;
 legend(text_metric);
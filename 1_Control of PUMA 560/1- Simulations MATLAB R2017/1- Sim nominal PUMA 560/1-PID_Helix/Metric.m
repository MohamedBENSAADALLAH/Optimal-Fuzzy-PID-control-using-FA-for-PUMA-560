function [ text_metric ] = Metric( error )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
RMSE=sqrt(mean((error).^2));
ISE=sum(error.*error);
str_RMSE =string(RMSE)
str_ISE =string(ISE)
text_metric = 'RMSE='+str_RMSE+' '+'ISE='+str_ISE
end


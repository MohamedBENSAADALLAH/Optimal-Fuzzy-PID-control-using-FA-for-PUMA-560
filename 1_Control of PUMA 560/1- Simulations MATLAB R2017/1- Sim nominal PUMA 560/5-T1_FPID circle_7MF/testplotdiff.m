diff =[]
for i =2:length(tout)
    diff=[diff;tout(i)-tout(i-1)];
end
plot(diff)
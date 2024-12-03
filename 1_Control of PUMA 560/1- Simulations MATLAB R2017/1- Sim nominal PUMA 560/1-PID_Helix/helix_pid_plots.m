close all; clc; clf ; 

%% Convert color code to 1-by-3 RGB array (0~1 each)
str1 = '#0072BD';
color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;
str2 = '#EDB120';
color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;
str3 = '#FF0000';
color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;

%% Plot Joints (j1,j2,j3) signals separetly
figure(1)
plot(j1_desired_helix_pid,'Color',color1,'LineWidth',1.5,...
    'MarkerSize',1); hold on; 
plot(j1_nominal_helix_pid,'Color',color2,'LineWidth',1.5,...
    'MarkerSize',5); hold on; 
plot (j1_perturbed_helix_pid,'--','Color',color3, 'LineWidth',1.5,...
    'MarkerSize',10); 
title('PID controller: 1st joint plots')
grid on 
xlabel('Time (s)')
ylabel('RMSE of J1')
legend('Desired J1','Nominal J1','Perturbed J1')

figure(2)
plot(j2_desired_helix_pid,'Color',color1,'LineWidth',1.5,...
    'MarkerSize',5); hold on; 
plot(j2_nominal_helix_pid,'Color',color2,'LineWidth',1.5,...
    'MarkerSize',10); hold on; 
plot (j2_perturbed_helix_pid,'--','Color',color3,'LineWidth',1.5,...
    'MarkerSize',5); 
title('PID controller: 2nd joint plots')
grid on
xlabel('Time (s)')
ylabel('RMSE of J2')
legend('Desired J2','Nominal J2','Perturbed J2')

figure(3)
plot(j3_desired_helix_pid,'Color',color1,'LineWidth',1.5,...
    'MarkerSize',5); hold on; 
plot(j3_nominal_helix_pid,'Color',color2,'LineWidth',1.5,...
    'MarkerSize',10); hold on; 
plot (j3_perturbed_helix_pid,'--','Color',color3,'LineWidth',1.5,...
    'MarkerSize',5); 
title('PID controller: 3rd joint plots')
grid on
xlabel('Time (s)')
ylabel('RMSE of J3')
legend('Desired J3','Nominal J3','Perturbed J3')

%% Plot xyz(3D) trajectories
figure(4)
%plot(xyz_des_helix_pid.Time,xyz_des_helix_pid.Data(:,1),xyz_des_helix_pid.Data(:,2),xyz_des_helix_pid.Data(:,3))
plot3(xyz_des_helix_pid.Data(:,1),xyz_des_helix_pid.Data(:,2),xyz_des_helix_pid.Data(:,3),'Color',color1,'LineWidth',1.5); hold on;
xlim([-0.2 0.8]); ylim([-0.6 0.2]); zlim([-0.6 0.6]) ; 

plot3(xyz_nom_helix_pid.Data(:,1),xyz_nom_helix_pid.Data(:,2),xyz_nom_helix_pid.Data(:,3),'Color',color2,'LineWidth',1.5); hold on ;
xlim([-0.2 0.8]); ylim([-0.6 0.2]); zlim([-0.6 0.6]) ; 

plot3(xyz_pert_helix_pid.Data(:,1),xyz_pert_helix_pid.Data(:,2),xyz_pert_helix_pid.Data(:,3),'--','Color',color3,'LineWidth',1.5);
xlim([-0.2 0.8]); ylim([-0.6 0.2]); zlim([-0.6 0.6]) ; 

title('PID controller: XYZ Helix Trajectory')
grid on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
legend('Desired trajectory','Nominal trajectory','Perturbed trajectory')



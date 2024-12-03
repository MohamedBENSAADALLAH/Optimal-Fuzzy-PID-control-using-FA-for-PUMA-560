function out3=TFC_T3(input3)

Ta=evalin('base','Ta');
warning off;

%input3=[0.5 0.6]

error=input3(1);
error_dot=input3(2);
tfc=newfis('tfc','sugeno');        % Sugeno Type FLC 


tfc = addvar(tfc,'input','error1',[-1 1]);
tfc = addvar(tfc,'input','error_dot1',[-1 1]);
tfc = addvar(tfc,'output','u1',[0 1]);
% MF for input1
tfc = addmf(tfc,'input',1,'NB','trapmf',[-10.5 -10 -1 -0.5]);         %Five Triangular Membership functions
tfc = addmf(tfc,'input',1,'NS','trimf',[-1 -0.5 0]);
tfc = addmf(tfc,'input',1,'ZO','trimf',[-0.5 0 0.5]);
tfc = addmf(tfc,'input',1,'PS','trimf',[0 0.5 1]);
tfc = addmf(tfc,'input',1,'PB','trapmf',[0.5 1 10 10.5]);
% MF for input2
tfc = addmf(tfc,'input',2,'NB','trapmf',[-10.5 -10 -1 -0.5]);         %Five Triangular Membership functions
tfc = addmf(tfc,'input',2,'NS','trimf',[-1 -0.5 0]);
tfc = addmf(tfc,'input',2,'ZO','trimf',[-0.5 0 0.5]);
tfc = addmf(tfc,'input',2,'PS','trimf',[0 0.5 1]);
tfc = addmf(tfc,'input',2,'PB','trapmf',[0.5 1 10 10.5]);
% MF for output1 
tfc = addmf(tfc,'output',1,'NB','constant',[-140]);                          %Five Triangular Membership functions
tfc = addmf(tfc,'output',1,'NS','constant',[-70]);
tfc = addmf(tfc,'output',1,'ZO','constant',[0]);
tfc = addmf(tfc,'output',1,'PS','constant',[70]);
tfc = addmf(tfc,'output',1,'PB','constant',[140]);


rule_matrix =[   1     1     1     1     1;
                 1     2     1     1     1;
                 1     3     2     1     1;
                 1     4     2     1     1;
                 1     5     3     1     1;
                 2     1     1     1     1;
                 2     2     1     1     1;
                 2     3     2     1     1;
                 2     4     3     1     1;
                 2     5     4     1     1;
                 3     1     2     1     1;
                 3     2     2     1     1;
                 3     3     3     1     1;
                 3     4     4     1     1;
                 3     5     4     1     1;
                 4     1     2     1     1;
                 4     2     3     1     1;
                 4     3     4     1     1;
                 4     4     5     1     1;
                 4     5     5     1     1;
                 5     1     3     1     1;
                 5     2     4     1     1;
                 5     3     4     1     1;
                 5     4     5     1     1;
                 5     5     5     1     1;
];

tfc=addrule(tfc,rule_matrix);
out3= evalfis([error error_dot],tfc);


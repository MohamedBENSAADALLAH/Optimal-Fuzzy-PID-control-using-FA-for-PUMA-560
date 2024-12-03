%________________________________________________________________________%
% Optimal PID-type Fuzzy Logic Controller source codes version 1.0      %
%                                                                        %
% Developed in MATLAB R2014b                                            %
%                                                                        %
% Author and programmer: Amin Noshadi                                   %
%                                                                        %
%         E-mail: Amin.Noshadi@live.vu.edu.au                            %
%                 Amin.Noshadi@gmail.com                                 %
%                                                                        %
%                                                                        %
% Please cite the Original Paper:                                      %
% A. Noshadi, J. Shi, W. S. Lee, P. Shi, and A. Kalam                    %
% Optimal PID-type Fuzzy Logic Controller for a Multi-Input Multi-Output %
% Active Magnetic Bearing System, Neural Computing and Applications,     %
% In Press, DOI: 10.1007/s00521-015-1996-7                                %
%                                                                        %
%________________________________________________________________________%

% Open the Simulink file "system2_2T.mdl" and run the m file "Fuzzy_Optimize_PSOMATLAB"

% This package is only for optimization of the scaling factors.
% It is assumed that the 5 triangular membership functions equally distributed with 25 rules.
% This may not be the best way to construct the Fuzzy Logic Controller (FLC), 
% because the FLC is reconstructed at every simulation run and hence it may be much slower than a fixed FLC structure.
% However, here we are able to optimize not only the input-output scaling factors,
% but also the the input/output membership function distribution and the rule-base. 
% Three sets of files will be uploaded eventually:
% 1- Only optimizes the scaling factors (T)
% 2- Only optimizes the distribution of the membership functions (S), 
% assuming that “T”s are already optimized.
% 3- Only optimizes the rule-base (R), assuming that “T”s and “S”s are already optimized.

% The number of membership functions, type of membership functions, 
% and the number of rule-bases can be modified in the TFC_T file.

% The main file “Fuzzy_Optimize_PSOMATLAB” uses particle swarm optimization (PSO) in MATLAB,
% however, other optimization algorithms will be added in the future. 

% lb and ub define the lower and upper bound on the search domain.


% The fitness function is defined as “fitness” and the description is given in
% the original paper “A. Noshadi, J. Shi, W. S. Lee, P. Shi, and A. Kalam,
% Optimal PID-type Fuzzy Logic Controller for a Multi-Input Multi-Output Active Magnetic Bearing System,
% Neural Computing and Applications, In Press, DOI: 10.1007/s00521-015-1996-7”. 

% Other objective functions will result in different solutions. 
% Furthermore, changing the values of Q and R in the objective function will also result in different solution.
% Try (Q=0.8,R=0.2) and (Q=0.2,R=0.8) to see the difference.


% It should be noted that the two input scaling factors (Ta(1)) and (Ta(2))
% are located in the “TFC_T” and the output scaling factors are located in the Simulink File. 


% Any question, don’t hesitate to contact me: Amin.Noshadi@live.vu.edu.au







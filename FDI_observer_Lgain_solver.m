% Copyright (c) 2026 [Wenbo Li]. All rights reserved.
%
% This code is the proprietary and confidential information of the author.
% It is provided for academic or evaluation purposes only.
% Unauthorized copying, modification, distribution, or commercial use of
% this code, in whole or in part, is strictly prohibited without prior
% written permission from the author.
%
% For permission to use this code, please contact: [wenbojlu@126.com]

% Solve for the observer gain matrix!!!

clc;
clear;
close all;

% Parameters of the system
dt  = 0.05; % Sampling time
tau = 0.6; % motor time constant

% System matrix-(equation (2))
A = [1,dt,0;
     0,1,dt;
     0,0,1-dt/tau];
B = [0;
     0;
     dt/tau];
C =  eye(3);

% Augmented system-(equation (14))
A1 = [A,B;zeros(1,3),1];
B1 = [B;0];
C1 = [C,zeros(3,1)];
W1 = [0 ;0 ;0 ;1];

 
%LMI-paramter-(equation (18))
setlmis([]);
F  = lmivar(1, [4 1]);    % psi
Fv = lmivar(2, [3 4]);    % lambda
gamma = lmivar(2, [1 1]);  

%LMI-constraints-(equation (18)) 
lmiterm([1 1 1 F], -1, 1);
lmiterm([1 1 1 0], eye(4));
lmiterm([1 1 3 F], A1', 1);
lmiterm([1 1 3 Fv],-C1', 1);      
lmiterm([1 2 2 gamma], -1, 1);
lmiterm([1 2 3 F],W1',1);
lmiterm([1 3 3 F], -1, 1);

lmiterm([-2 1 1 F], 1, 1);  % psi > 0

% 获取LMI系统
LMIs = getlmis;

[tmin, xfeas] = feasp(LMIs);

if  tmin < 0
    F_sol  = dec2mat(LMIs, xfeas, F);     
    Fv_sol = dec2mat(LMIs, xfeas, Fv);  
    disp('It has a feasible solution');
else
    disp('It has not a feasible solution');
end

L = F_sol^-1*Fv_sol';
 



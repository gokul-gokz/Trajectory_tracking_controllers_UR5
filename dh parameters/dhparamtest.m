%% A code to check that the DH-Parameters are correct
% Using the Robotix Toolbox for Matlab (included in the folder)

close all
clear
clc

L1 = Link('d', 0.08916, 'a', 0, 'alpha', pi/2);
L2 = Link('d', 0, 'a', -0.425, 'alpha', 0);
L3 = Link('d', 0, 'a', -0.39225, 'alpha', 0);
L4 = Link('d', 0.10915, 'a', 0, 'alpha', pi/2);
L5 = Link('d', 0.09456, 'a', 0, 'alpha', -pi/2);
% tool frame
L6 = Link('d', 0.0823, 'a', 0, 'alpha', 0);

bot = SerialLink([L1 L2 L3 L4 L5 L6], 'name', 'ur5');

% draw the robot (do not change the values here)
bot.plot([0 0 0 0 0 0]);

%% compute the FK

syms q1 q2 q3 q4 q5 q6

dhparams = [   q1,   0.08916,        0,             pi/2; ...
               q2,   0,             -0.425,         0; ...
               q3,   0,             -0.39225,       0; ...
               q4,   0.10915,        0,             pi/2; ...
               q5,   0.09456,        0,             -pi/2; ...
               q6,   0.0823,         0,             0];

% compute the forward kinematics using the DH-parameters
T = dhparam2matrices(dhparams(:,1), dhparams(:,2), dhparams(:,3), dhparams(:,4), false);

T_of_interest = T(:,:,6);
FK = T_of_interest(1:3, 4);

FK_f = matlabFunction(FK, 'File','computeFK','Vars',{[q1; q2; q3; q4; q5; q6]});
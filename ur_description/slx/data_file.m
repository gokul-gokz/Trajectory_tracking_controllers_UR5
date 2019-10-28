%% Initialization
% must be inside the ur_description/slx folder to run
addpath('../meshes/ur5/collision');
addpath('../meshes/ur5/visual');
addpath('../../dynamic modeling');
addpath('../../dh parameters');

%% Damping
% specifies the damping in 
data.damping_coefficient = 0; % in N*m/(rad/s)


% generate the trajectories
[trajectories, parameters] = generatePolynomialManipulatorTrajectory(5, [0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0], [0, pi/2, 0, 0, 0, 0; 0 0 0 0 0 0;0 0 0 0 0 0], 10);

% fetch the parameters for each joint
a1 = parameters(:,1);
a2 = parameters(:,2);
a3 = parameters(:,3);
a4 = parameters(:,4);
a5 = parameters(:,5);
a6 = parameters(:,6);

T = [0:0.01:10];
q_d = vpa(subs(trajectories(1,2),T),5);

plot(T, q_d);

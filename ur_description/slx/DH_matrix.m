function [T]= DH_matrix(a, alpha, theta, d)

Ta=[1 0 0 a; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Rz= [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
Rx=[1 0 0 0; 0 cos(alpha) -sin(alpha) 0; 0 sin(alpha) cos(alpha) 0; 0 0 0 1];
Td=[1 0 0 0; 0 1 0 0; 0 0 1 d; 0 0 0 1];
T= Rz * Td * Ta * Rx;   
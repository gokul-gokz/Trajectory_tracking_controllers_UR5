clc; clear all;
syms q1 q2 q3 q4 q5 q6
g = 9.8;
m1 = 3.7;
m2 = 8.393;
m3 = 2.33;
m4 = 1.219;
m5 = 1.219;
m6 = 0.1879;
m = [m1, m2, m3, m4, m5, m6];
%% Link lengths
l1 = 130*1e-3;
l2 = 542.8*1e-3;
l3 = 489.3*1e-3;
l4 = 105*1e-3;
l5 = 105*1e-3;
l6 = 36*1e-3;
l = [l1, l2, l3, l4, l5, l6];

%% Moment of inertia
I1 = (1/3)*m1*l1^2;
I2 = (1/3)*m2*l2^2;
I3 = (1/3)*m3*l3^2;
I4 = (1/3)*m4*l4^2;
I5 = (1/3)*m5*l5^2;
I6 = (1/3)*m6*l6^2;
I = [I1, I2, I3, I4, I5, I6];

%% Center of mass
r1 = [0;-0.02561;0.00193];
r2 = [0.2125; 0; 0.11336];
r3 = [0.15; 0; 0.0265];
r4 = [0; -0.0018; 0.01634];
r5 = [0; 0.0018; 0.01634];
r6 = [0; 0; -0.001159];
r = [r1, r2, r3, r4, r5, r6];

%% Transformation Matrices Using DH Parameters (a, alpha, theta, d)

T10 = DH_matrix(0, pi/2, q1, 0.08916);
T21 = DH_matrix(-0.425, 0, q2, 0);
T32 = DH_matrix(-0.39225, 0, q3, 0);
T43 = DH_matrix(0, pi/2, q4, 0.10915);
T54 = DH_matrix(0, -pi/2, q5, 0.09456);
T65 = DH_matrix(0, 0, q6, 0.0823);

%% Transformation Matrices wrt Base frames
T20 = T10 * T21;
T30 = T20 * T32;
T40 = T30 * T43;
T50 = T40 * T54;
T60 = T50 * T65;


%% Jacobians
% Jv
%% Link 1
J1 = diff(T10*[0;0;0;1], q1);
J1 = J1(1:3);
J2 = diff(T10*[0;0;0;1], q2);
J2 = J2(1:3);
J3 = diff(T10*[0;0;0;1], q3);
J3 = J3(1:3);
J4 = diff(T10*[0;0;0;1], q4);
J4 = J4(1:3);
J5 = diff(T10*[0;0;0;1], q5);
J5 = J5(1:3);
J6 = diff(T10*[0;0;0;1], q6);
J6 = J6(1:3);
Jv1 = vpa([J1, J2, J3, J4, J5, J6])
Jv(1,:,:) = Jv1;
%% Link 2
J1 = diff(T20*[0;0;0;1], q1);
J1 = J1(1:3);
J2 = diff(T20*[0;0;0;1], q2);
J2 = J2(1:3);
J3 = diff(T20*[0;0;0;1], q3);
J3 = J3(1:3);
J4 = diff(T20*[0;0;0;1], q4);
J4 = J4(1:3);
J5 = diff(T20*[0;0;0;1], q5);
J5 = J5(1:3);
J6 = diff(T20*[0;0;0;1], q6);
J6 = J6(1:3);
Jv2 = vpa([J1, J2, J3, J4, J5, J6]);
Jv(2,:,:) = Jv2;
%% Link 3
J1 = diff(T30*[0;0;0;1], q1);
J1 = J1(1:3);
J2 = diff(T30*[0;0;0;1], q2);
J2 = J2(1:3);
J3 = diff(T30*[0;0;0;1], q3);
J3 = J3(1:3);
J4 = diff(T30*[0;0;0;1], q4);
J4 = J4(1:3);
J5 = diff(T30*[0;0;0;1], q5);
J5 = J5(1:3);
J6 = diff(T30*[0;0;0;1], q6);
J6 = J6(1:3);
Jv3 = vpa([J1, J2, J3, J4, J5, J6]);
Jv(3,:,:) = Jv3;
%% Link 4
J1 = diff(T40*[0;0;0;1], q1);
J1 = J1(1:3);
J2 = diff(T40*[0;0;0;1], q2);
J2 = J2(1:3);
J3 = diff(T40*[0;0;0;1], q3);
J3 = J3(1:3);
J4 = diff(T40*[0;0;0;1], q4);
J4 = J4(1:3);
J5 = diff(T40*[0;0;0;1], q5);
J5 = J5(1:3);
J6 = diff(T40*[0;0;0;1], q6);
J6 = J6(1:3);
Jv4 = vpa([J1, J2, J3, J4, J5, J6])
Jv(4,:,:) = Jv4;
%% Link 5
J1 = diff(T50*[0;0;0;1], q1);
J1 = J1(1:3);
J2 = diff(T50*[0;0;0;1], q2);
J2 = J2(1:3);
J3 = diff(T50*[0;0;0;1], q3);
J3 = J3(1:3);
J4 = diff(T50*[0;0;0;1], q4);
J4 = J4(1:3);
J5 = diff(T50*[0;0;0;1], q5);
J5 = J5(1:3);
J6 = diff(T50*[0;0;0;1], q6);
J6 = J6(1:3);
Jv5 = vpa([J1, J2, J3, J4, J5, J6]);
Jv(5,:,:) = Jv5;
%% Link 6
J1 = diff(T60*[0;0;0;1], q1);
J1 = J1(1:3);
J2 = diff(T60*[0;0;0;1], q2);
J2 = J2(1:3);
J3 = diff(T60*[0;0;0;1], q3);
J3 = J3(1:3);
J4 = diff(T60*[0;0;0;1], q4);
J4 = J4(1:3);
J5 = diff(T60*[0;0;0;1], q5);
J5 = J5(1:3);
J6 = diff(T60*[0;0;0;1], q6);
J6 = J6(1:3);
Jv6 = vpa([J1, J2, J3, J4, J5, J6]);
Jv(6,:,:) = Jv6;
% Jw
Jw1 = T10*[0;0;1;0];
Jw1 = Jw1(1:3);
Jw2 = T20*[0;0;1;0];
Jw2 = Jw2(1:3);
Jw3 = T30*[0;0;1;0];
Jw3 = Jw3(1:3);
Jw4 = T40*[0;0;1;0];
Jw4 = Jw4(1:3);
Jw5 = T50*[0;0;1;0];
Jw5 = Jw2(1:3);
Jw6 = T60*[0;0;1;0];
Jw6 = Jw6(1:3);
Jw = vpa([Jw1, Jw2, Jw3, Jw4, Jw5, Jw6])
% J = [Jv;Jw];
R(1,:,:) = T10(1:3,1:3);
R(2,:,:) = T20(1:3,1:3);
R(3,:,:) = T30(1:3,1:3);
R(4,:,:) = T40(1:3,1:3);
R(5,:,:) = T50(1:3,1:3);
R(6,:,:) = T60(1:3,1:3);
%% Dynamics
M = zeros(6,6);
for i=1:6
    
    R_ = reshape(R(i,:,:),[3,3])
    J = reshape(Jv(i,:,:), [3, 6])
    transpose(J)*J
    size(Jw(:,i))
    transpose(Jw(:,i))*R_*I(i)*transpose(R_)*Jw(:,i)
    M = vpa(M + m(i)*transpose(J)*J +transpose(Jw(:,i))*Rg _*I(i)*transpose(R_)*Jw(:,i));
end
M
%% Finding Transformation Matrix from DH parameters
function [T]= DH_matrix(a, alpha, theta, d)
    Ta=[1 0 0 a; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    Rz= [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
    Rx=[1 0 0 0; 0 cos(alpha) -sin(alpha) 0; 0 sin(alpha) cos(alpha) 0; 0 0 0 1];
    Td=[1 0 0 0; 0 1 0 0; 0 0 1 d; 0 0 0 1];
    T= Rz * Td * Ta * Rx;  
end

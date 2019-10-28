clc; clear all;
syms q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) t1 t2 t3 t1_d t2_d t3_d t1_dd t2_dd t3_dd t g m1 m2 m3;
%      l1 l2 l3 I1 I2 I3
% m1 = 3.7;
% m2 = 8.393;
% m3 = 2.33;
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
l = sym([l1, l2, l3, l4, l5, l6]);

%% Moment of inertia
I1 = (1/12)*m1*l1^2;
I2 = (1/12)*m2*l2^2;
I3 = (1/12)*m3*l3^2;
I4 = (1/12)*m4*l4^2;
I5 = (1/12)*m5*l5^2;
I6 = (1/12)*m6*l6^2;
I = [I1, I2, I3, I4, I5, I6];

%% Center of mass
r1 = [0;-0.02561;0.00193];
r2 = [0.2125; 0; 0.11336];
r3 = [0.15; 0; 0.0265];
r4 = [0; -0.0018; 0.01634];
r5 = [0; 0.0018; 0.01634];
r6 = [0; 0; -0.001159];
r = sym([r1, r2, r3, r4, r5, r6]);

%% Transformation Matrices Using DH Parameters (a, alpha, theta, d)

T10 = DH_matrix(0, sym( pi/2), q1, sym(0.08916));
T21 = DH_matrix(sym(-0.425), 0, q2, 0);
T32 = DH_matrix(sym(-0.39225), 0, q3, 0);
T43 = DH_matrix(0, sym(pi/2), 0, sym(0.10915));
T54 = DH_matrix(0, sym(-pi/2), 0, sym(0.09456));
T65 = DH_matrix(0, 0, 0, sym(0.0823));

%% Transformation Matrices wrt Base frames
T20 = T10 * T21;
T30 = T20 * T32;
T40 = T30 * T43;
T50 = T40 * T54;
T60 = T50 * T65;

%% Positions
x1 = T10*[0;0;0;1];
x2 = T20*[0;0;0;1];
x3 = T30*[0;0;0;1];
x4 = T40*[0;0;0;1];
x5 = T50*[0;0;0;1];
x6 = T60*[0;0;0;1];

%% Velocities
v1 = diff(x1,t);
v2 = diff(x2,t);
v3 = diff(x3,t);
v4 = diff(x4,t);
v5 = diff(x5,t);
v6 = diff(x6,t);

%% Kinetic Energy
K1 = 0.5*m(1)*transpose(v1)*v1 + 0.5*I(1)*diff(q1,t)^2;
K2 = 0.5*m(2)*transpose(v2)*v2 + 0.5*I(2)*diff(q2,t)^2;
K3 = 0.5*m(3)*transpose(v3)*v3 + 0.5*I(3)*diff(q3,t)^2;
K4 = 0.5*m(4)*transpose(v4)*v4 + 0.5*I(4)*diff(q4,t)^2;
K5 = 0.5*m(5)*transpose(v5)*v5 + 0.5*I(5)*diff(q5,t)^2;
K6 = 0.5*m(6)*transpose(v6)*v6 + 0.5*I(6)*diff(q6,t)^2;
K = K1+K2+K3+K4+K5+K6;
%% Potential Energy
P1 = m(1)*g*[0 0 1 0]*x1;
P2 = m(2)*g*[0 0 1 0]*x2;
P3 = m(3)*g*[0 0 1 0]*x3;
P4 = m(4)*g*[0 0 1 0]*x4;
P5 = m(5)*g*[0 0 1 0]*x5;
P6 = m(6)*g*[0 0 1 0]*x6;
P = P1+P2+P3+P4+P5+P6;


%% Lagrange
L = K-P;
exp1 = subs(L, [q1(t), q2(t), q3(t), diff(q1(t),t), diff(q2(t),t),diff(q3(t),t)],[t1,t2,t3,t1_d,t2_d,t3_d]);
exp2 = [diff(exp1, t1); diff(exp1, t2); diff(exp1, t3)];
exp3 = [diff(exp1, t1_d);diff(exp1, t2_d);diff(exp1, t3_d)];
exp4 = diff(subs(exp3,[t1,t2,t3,t1_d,t2_d,t3_d], [q1(t), q2(t), q3(t), diff(q1(t),t), diff(q2(t),t),diff(q3(t),t)]),t);
exp = subs(exp4-exp1,[q1(t), q2(t), q3(t), diff(q1(t),t), diff(q2(t),t),diff(q3(t),t),  diff(q1,t,t), diff(q2,t,t), diff(q3,t,t)],...
      [t1,t2,t3,t1_d,t2_d,t3_d,t1_dd,t2_dd,t3_dd]);
  
%% Torque
Torque = simplify(expand(exp));
th = [t1, t2, t3];
dth = [t1_d, t2_d, t3_d];

%% Mass, Coriolis Matrices
M = vpa(equationsToMatrix(transpose(Torque), [t1_dd, t2_dd, t3_dd]))
C = sym(zeros(3,3));
for i=1:3
    for j=1:3
        C(i,j) = 0;
        C(i,j)= (1/2)*(diff(M(i,j),th(1))+diff(M(i,1),th(j))-diff(M(1,j),th(i))*dth(1)+...
                            (diff(M(i,j),th(2))+diff(M(i,2),th(j))-diff(M(2,j),th(i)))*dth(2)+...
                            (diff(M(i,j),th(3))+diff(M(i,3),th(j))-diff(M(3,j),th(i)))*dth(3));
    end
end
G = vpa(equationsToMatrix(transpose(Torque), g))
C = vpa(C)

%% Jacobian
%% Jacobian
posx = subs([1 0 0 0]*x6, [q1, q2, q3], [t1,t2,t3]);
posy = subs([0 1 0 0]*x6, [q1, q2, q3], [t1,t2,t3]); 
posz = subs([0 0 1 0]*x6, [q1, q2, q3], [t1,t2,t3]);
J = vpa([diff(posx, t1), diff(posx, t2), diff(posx, t3);...
     diff(posy, t1), diff(posy, t2), diff(posy, t3);...
     diff(posz, t1), diff(posz, t2), diff(posz, t3)])
dJ = vpa(subs(diff(subs(J, [t1,t2,t3],[q1, q2, q3]), t), [q1(t), q2(t), q3(t), diff(q1(t),t), diff(q2(t),t),diff(q3(t),t)],[t1,t2,t3,t1_d,t2_d,t3_d]))


%% Finding Transformation Matrix from DH parameters
function [T]= DH_matrix(a, alpha, theta, d)
    Ta=[1 0 0 a; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    Rz= [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
    Rx=[1 0 0 0; 0 cos(alpha) -sin(alpha) 0; 0 sin(alpha) cos(alpha) 0; 0 0 0 1];
    Td=[1 0 0 0; 0 1 0 0; 0 0 1 d; 0 0 0 1];
    T= Rz * Td * Ta * Rx;  
end

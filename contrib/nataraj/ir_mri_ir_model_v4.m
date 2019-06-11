% script ir_mri_IR_model_v4.m
% deriving a more thorough Inversion-Recovery Spin-Echo (IR SE) model
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%	1.1	2016-02-12	original
%	1.2	2016-05-27	added gradient; generate signal function
%	1.3	2016-06-22	testing whether alternating rf excitation axis changes anything
%	1.4	2016-06-26	accounting for imperfect inversion
% 2016-08-10 JF changed name and white space :)

% variable definitions
syms M0 T1 T2 inveff kap TR TI TE t positive;
syms a flip_inv flip_ex flip_ref dw real;

% assumptions
assume((TR > TI) & (TI > TE) & (TE > 0) & (T1 > T2));

% add spatial variation to flip angles
a_inv = inveff * flip_inv;
a_ex = kap * flip_ex;
a_ref = kap * flip_ref;


%% define matrix operations
% excitation about x'
Rx(a) =     [1,             0,              0;...
             0,             cos(a),         sin(a);...
             0,             -sin(a),        cos(a)];

% excitation about y'
Ry(a) =     [cos(a),        0,              sin(a);...
             0,             1,              0;...
             -sin(a),       0,              cos(a)];

% precession and relaxation
E(t) =      [exp(-t/T2),    0,              0;...
             0,             exp(-t/T2),     0;...
             0,             0,              exp(-t/T1)];
Rz(t) =     [cos(dw*t),     sin(dw*t),      0;...
             -sin(dw*t),    cos(dw*t),      0;...
             0,             0,              1];
A(t) = Rz(t) * E(t);

% ideal spoiling
S =         [0,             0,              0;...
             0,             0,              0;...
             0,             0,              1];

% mz recovery
b(t) =      [0;             0;              M0*(1-exp(-t/T1))];


%% matrix equation solution
syms m1_TI_TE_x m1_TI_TE_y m1_TI_TE_z;
m1_TI_TE = [m1_TI_TE_x; m1_TI_TE_y; m1_TI_TE_z];

m2_z_m = A(TR-(TI+TE))*S*m1_TI_TE + b(TR-(TI+TE));
m2_z_p = S*Rx(a_inv)*m2_z_m; % Rx or Ry doesn't matter
m2_TI_m = A(TI)*m2_z_p + b(TI);
% m2_TI_p = Rx(+a_ex)*m2_TI_m; % m1 is excited by -a_ex
m2_TI_p = Rx(-a_ex)*m2_TI_m; % m1 is excited by +a_ex
m2_TI_hTE_m = A(TE/2)*m2_TI_p + b(TE/2);
m2_TI_hTE_p = Ry(a_ref)*m2_TI_hTE_m;
m2_TI_TE = A(TE/2)*m2_TI_hTE_p + b(TE/2);

m3_z_m = A(TR-(TI+TE))*S*m2_TI_TE + b(TR-(TI+TE));
m3_z_p = S*Rx(a_inv)*m3_z_m; % Rx or Ry doesn't matter
m3_TI_m = A(TI)*m3_z_p + b(TI);
% m3_TI_p = Rx(-a_ex)*m3_TI_m; % m1 is excited by -a_ex
m3_TI_p = Rx(+a_ex)*m3_TI_m; % m1 is excited by +a_ex
m3_TI_hTE_m = A(TE/2)*m3_TI_p + b(TE/2);
m3_TI_hTE_p = Ry(a_ref)*m3_TI_hTE_m;
m3_TI_TE = A(TE/2)*m3_TI_hTE_p + b(TE/2);

[B1, y1] = equationsToMatrix(m3_TI_TE == m1_TI_TE, m1_TI_TE);
m1_sol = simplify(linsolve(B1, y1));
m2_sol = simplify(subs(m2_TI_TE,...
	[m1_TI_TE_x, m1_TI_TE_y, m1_TI_TE_z],...
	[m1_sol(1),  m1_sol(2),  m1_sol(3)]));

% display ideal case
pretty(simplify(subs(m1_sol, [kap, flip_ex, flip_ref], [1, pi/2, pi])));
pretty(simplify(subs(m2_sol, [kap, flip_ex, flip_ref], [1, pi/2, pi])));

% received signals
s1_TI_TE_xy = m1_sol(1) + 1i*m1_sol(2);
s2_TI_TE_xy = m2_sol(1) + 1i*m2_sol(2);

% generate signal functions
if ~exist('IR_gen_v3.m', 'file')
	tmp = matlabFunction(s1_TI_TE_xy, s2_TI_TE_xy,...
	'file', 'IR_gen_v3.m',...
	'vars', [M0 T1 T2 TR TI TE flip_inv flip_ex flip_ref inveff kap dw]);
end

% row gradient of received signals w.r.t. independent variables, x
s1_TI_TE_xy_gradx = simplify(jacobian(s1_TI_TE_xy, [M0 T1 T2 inveff])); % [1 L]
s2_TI_TE_xy_gradx = simplify(jacobian(s2_TI_TE_xy, [M0 T1 T2 inveff])); % [1 L]

% generate gradient function
if ~exist('IR_gradx_gen_v2.m', 'file')
	tmp = matlabFunction(s1_TI_TE_xy_gradx, s2_TI_TE_xy_gradx,...
	'file', 'IR_gradx_gen_v2.m',...
	'vars', [M0 T1 T2 TR TI TE flip_inv flip_ex flip_ref inveff kap dw]);
end

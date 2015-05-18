%
% new_dgm.m
%
% Copyright (C) 2015 Mathieu Gaborit (matael) <mathieu@matael.org>
%
%
% Distributed under WTFPL terms
%
%            DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%                    Version 2, December 2004
%
% Copyright (C) 2004 Sam Hocevar <sam@hocevar.net>
%
% Everyone is permitted to copy and distribute verbatim or modified
% copies of this license document, and changing it is allowed as long
% as the name is changed.
%
%            DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
%
%  0. You just DO WHAT THE FUCK YOU WANT TO.
%

clear all;
close all;


f = 500;
l = 1;

rho = 1.241;
c = 343;
Z0  = rho*c;

omega = 2*pi*f;
k = omega/c;

Ax = [
	0, 1/rho;
	rho*c^2, 0
];

% diagonalization
P = [
	1, 1;
	Z0, -Z0
];
Pp = P(:,1);
Pm = P(:,2);

L = [c, 0; 0, -c];
Lm = -c;
Lp = c;


Q = [
	0.5, 1/(2*Z0);
	0.5, -1/(2*Z0)
];
Qp = Q(1,:);
Qm = Q(2,:);


Uu = @(x) [
	exp(-j*k*x)/Z0, -exp(j*k*x)/Z0;
	exp(-j*k*x), exp(j*k*x);
];

Uv = @(x) [
	exp(-j*k*x)*Z0, -exp(j*k*x)*Z0;
	exp(-j*k*x), exp(j*k*x)
];

% Boundary conditions
Rm = -1; % right BC

LH_mat = Uv(l)*P*L*Q*(Pp+Pm*Rm)*Qp*Uu(l);
LH_mat = LH_mat - Uv(0)*P*L*Q*Pm*Qm*Uu(0);

RH_vec = Uv(0)*P*L*Q*Pp;


U_vec = LH_mat\RH_vec;

x_v = 0:0.005:l;
pv = zeros(2,length(x_v));
for idx=1:length(x_v)
	pv(:,idx) = Uu(x_v(idx))*U_vec;
end

p = pv(2,:)/Z0;



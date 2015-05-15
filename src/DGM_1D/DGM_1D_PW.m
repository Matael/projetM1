%
% DGM_1D_PW.m
%
% Copyright (C) 2015 Mathieu Gaborit (matael) <mathieu@matael.org>
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

clear all;
close all;

% frequency
f = 500; % Hz

% medium
rho = 1.241;
c= 343;
Z0 = rho*c;

% geometry
l = 1;

% helpers
omega = 2*pi*f;
k = omega/c;

% Flux matrix
F = [0 , 1/rho ; rho*c^2 , 0];

L = [c, 0 ; 0, -c];
P = [1, 1; Z0, -Z0];
Q = inv(P);

% DEBUG
disp('Vérification de la diagonalisation : P*L*Q - F =')
P*L*Q-F


% Boundary condition matrices
Cm = [P(:,1) , [-1/Z0; 1]];
Cp = [[1/Z0; 1], P(:,2)];

disp('Calcul de R tilde : ')
Rtilde = inv(Cm)*Cp


val_pos = real(diag(L)) >= 0;

Lm = L;
Lm(:,find(val_pos)) = zeros(2,1);
Lp = L-Lm;
Qm = Q;
Qm(:,find(val_pos)) = zeros(2,1);
Qp = Q-Qm;

xm = l/2;
Uu = @(x) [
	exp(-j*k*x)/Z0 , -exp(j*k*x)/Z0;
	exp(-j*k*x) , exp(j*k*x)
];
% Uvt = @(x) transpose(conj(Uu(x)));
Uvt = @(x) transpose(Uu(x));

% matrix on the left-hand side
% LH_mat = -Uvt(l)*P*L*Qm*Uu(l)+Uvt(0)*P*(Lp*Rtilde(1,2)+Lm)*Qm*Uu(0);
% RH_vec = -Uvt(0)*P*Lp*Rtilde(1,1)*ones(2,1);
LH_mat = Uvt(l)*P*L*Qm*Uu(l)-Uvt(0)*P*(Lp*Rtilde(1,2)+Lm)*Qm*Uu(0);
RH_vec = Uvt(0)*P*Lp*Rtilde(1,1)*ones(2,1);

amplitudes = LH_mat\RH_vec;

R = Rtilde(2,1) + Rtilde(2,2)*Qm*Uu(0)*amplitudes(1:2);
angle(R)

Zi = Z0/(j*tan(k*l));
R_ana = angle((Zi-Z0)/(Zi+Z0))

figure;
x_v = (0:.005:l);
for idx=1:length(x_v)
	xi = x_v(idx);
	Uu_xi = Uu(xi);
	p(idx) = Uu_xi(2,:)*amplitudes(1:2);
end
plot(x_v, p, 'LineWidth', 2);

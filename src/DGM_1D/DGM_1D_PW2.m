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
F = [0 -1/rho ; rho*c^2 0];
S = [Z0 -Z0; 1 1]*[1;0];

% Boundary condition matrices
Cm = [1 Z0 ; -Z0 -1];
Cp = [Z0 -1 ; 1 -Z0];
Rtilde = -inv(Cm)*Cp;


[P, L] = eig(F);
Q = inv(P);

val_pos = imag(diag(L)) >= 0;
val_neg = ones(length(diag(L)),1)-val_pos;

Lm = L;
Lm(:,find(val_pos)) = zeros(2,1);
Lp = L-Lm;
Qp = Q;
Qp(:,find(val_pos)) = zeros(2,1);

Sm = inv(Cm)*S;

xm = l/2;
Uu = [
	exp(j*k*xm)/Z0 , -exp(-j*k*xm)/Z0;
	exp(j*k*xm) , exp(-j*k*xm)
];

Uv = Uu;

% matrix on the left-hand side
LH_mat = Uv'*P*(Lp-Lm*Rtilde)*Qp*Uu;
RH_vec = -Uv'*P*Lm*Sm;


amplitudes = LH_mat\RH_vec;

%
% 1D_FEM_PW.m
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

% fluid params
rho = 1.2141;
c = 343;
Zc = rho*c;

% frequency
f = 500; % Hz
k = 2*pi*f/c;

% geometrical
L = 1;

% Analytical solution
Zi = Zc/(j*tan(k*L));
R_ana = (Zi-Zc)/(Zi+Zc);

N_elem_max=400;

err = zeros(N_elem_max, 1);
R_vect_ref = zeros(N_elem_max, 1);
R_vect_test = zeros(N_elem_max, 1);
err_R_ref = zeros(N_elem_max, 1);
err_R_test = zeros(N_elem_max, 1);

for N_elem=1:N_elem_max
	% Discretization
	h = L/(2*N_elem);

	% elementary matrices
	M_elem = h/6*[2 1; 1 2];
	K_elem = 1/h*[1 -1; -1 1];

	% global matrices construction
	M_global = sparse(2*N_elem+1,2*N_elem+1);
	K_global = sparse(2*N_elem+1,2*N_elem+1);
	for i=1:N_elem
		M_global(i:i+1,i:i+1) = M_global(i:i+1,i:i+1) + M_elem;
		K_global(i:i+1,i:i+1) = K_global(i:i+1,i:i+1) + K_elem;
	end

	% global matrix
	A = sparse(2*N_elem+2,2*N_elem+2);
	A(1:2*N_elem+1,1:2*N_elem+1) = k^2*M_global - K_global;

	A(2*N_elem+2,1) = 1;
	A(2*N_elem+2,2*N_elem+2) = -1;
	A(1,2*N_elem+2) = -j*k;

	b = zeros(2*N_elem+2,1);
	b(1) = -j*k;
	b(2*N_elem+2) = 1;

	% resolution
	x = A\b;

	p_ref = x(1:2*N_elem+1);
	R_vect_ref(N_elem) = x(2*N_elem+2);

	% with DGM like formulation
	A = sparse(2*N_elem+1,2*N_elem+1);
	A = k^2*M_global - K_global;

	b = zeros(2*N_elem+1,1);
	b(1) = -j*k;
	A(1,1) = A(1,1) - 0.5*(j*k-1/h);
	A(1,2) = A(1,2) - 0.5/h;

	p_test = A\b;
	R_vect_test(N_elem) = 1/2*((1-1/(j*h*k))*p_test(1)+1/(j*h*k)*p_test(2));

	err(N_elem) = norm(p_ref-p_test, 2)/norm(p_ref,2);
	err_R_ref(N_elem) = norm(angle(R_vect_ref(N_elem))-angle(R_ana),2)/norm(angle(R_ana),2);
	err_R_test(N_elem) = norm(angle(R_vect_test(N_elem))-angle(R_ana),2)/norm(angle(R_ana),2);
end

figure;
plot(real(p_ref), '+r', 'LineWidth', 2);
hold on;
plot(real(p_test), 'ob', 'LineWidth', 2);

figure;
loglog(1:N_elem_max, err, '+r', 'LineWidth', 2)
xlabel('Nb DOF')
grid on;

figure;
loglog(1:N_elem_max, err_R_ref, 'b', 'LineWidth', 2)
hold on;
loglog(1:N_elem_max, err_R_test, 'r', 'LineWidth', 2)

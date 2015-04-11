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

N_elem_max = 100;
R_vect = zeros(N_elem_max,1);
for N_elem=1:100
	% Discretization
	h = L/N_elem;

	% elementary matrices
	M_elem = h/6*[2 1; 1 2];
	K_elem = 1/h*[1 -1; -1 1];

	% global matrices construction
	M_global = sparse(N_elem+1,N_elem+1);
	K_global = sparse(N_elem+1,N_elem+1);
	for i=1:N_elem
		M_global(i:i+1,i:i+1) = M_global(i:i+1,i:i+1) + M_elem;
		K_global(i:i+1,i:i+1) = K_global(i:i+1,i:i+1) + K_elem;
	end

	% global matrix
	A = sparse(N_elem+2,N_elem+2);
	A(1:N_elem+1,1:N_elem+1) = k^2*M_global - K_global;

	A(N_elem+2,1) = 1;
	A(N_elem+2,N_elem+2) = -1;
	A(1,N_elem+2) = -j*k;

	b = zeros(N_elem+2,1);
	b(1) = -j*k;
	b(N_elem+2) = 1;

	% resolution
	x = A\b;

	p = x(1:N_elem+1);
	R = x(N_elem+2);

	R_vect(N_elem) = R;
end

% Analytical solution
Zi = Zc/(j*tan(k*L));
R_ana = (Zi-Zc)/(Zi+Zc);

% error
err = zeros(N_elem_max, 1);
for i=1:N_elem_max
	err(i) = norm(angle(R_vect(i)) - angle(R_ana),2)/norm(angle(R_ana),2);
end


% figures
figure; % error (LogLog)
loglog(1:N_elem_max , err, 'LineWidth', 2)
xlabel('Number of elements')
ylabel('Relative Error')
grid on;
set(gca, 'xminorgrid', 'off');
set(gca, 'yminorgrid', 'off');

figure; % R FEM & R ana (angle)
plot(1:N_elem_max, angle(R_vect)+(angle(R_vect)<0)*2*pi, 'r+');
hold on;
plot([1 N_elem_max], [1 1]*angle(R_ana), 'k:', 'LineWidth', 2);
xlabel('Number of elements')
ylabel('Phase of R')
legend('FEM', 'Analytical')
grid on;

figure; % R in z-zplane
plot(real(R_vect), imag(R_vect), '+b')
hold on;
plot(real(R_ana), imag(R_ana), 'xor', 'LineWidth', 2)
grid on;



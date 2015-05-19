%
% FEM_1D_PW_quad.m
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

% load utils
addpath('../utils');
load_utils

% fluid params
rho = 1.2141;
c = 343;
Zc = rho*c;

% frequency
f = 1000; % Hz
k = 2*pi*f/c;

% geometrical
L = 1;

N_elem_max = 400;
R_vect = zeros(N_elem_max,1);
for N_elem=1:N_elem_max

	% Discretization
	h = L/N_elem;

	% form functions
	phi1 = @(x) ((h-2*x)*(h-x))/(h^2);
	phi2 = @(x) -4*x*(x-h)/h^2;
	phi3 = @(x) x*(2*x-h)/h^2;

	Iphi1_2 = intg_fun(mulfun(phi1, phi1), 0, h);
	Iphi2_2 = intg_fun(mulfun(phi2, phi2), 0, h);
	Iphi3_2 = Iphi1_2;
	Iphi1phi2 = intg_fun(mulfun(phi1, phi2), 0, h);
	Iphi1phi3 = intg_fun(mulfun(phi1, phi3), 0, h);
	Iphi2phi3 = intg_fun(mulfun(phi2, phi3), 0, h);

	drv_phi1 = @(x) (4*x-3*h)/(h^2);
	drv_phi2 = @(x) -4*(2*x-h)/h^2;
	drv_phi3 = @(x) (4*x-h)/h^2;

	Idrv_phi1_2 = intg_fun(mulfun(drv_phi1, drv_phi1), 0, h);
	Idrv_phi2_2 = intg_fun(mulfun(drv_phi2, drv_phi2), 0, h);
	Idrv_phi3_2 = Idrv_phi1_2;
	Idrv_phi1drv_phi2 = intg_fun(mulfun(drv_phi1, drv_phi2), 0, h);
	Idrv_phi1drv_phi3 = intg_fun(mulfun(drv_phi1, drv_phi3), 0, h);
	Idrv_phi2drv_phi3 = intg_fun(mulfun(drv_phi2, drv_phi3), 0, h);

	% elementary matrices
	M_elem = [
		Iphi1_2 , Iphi1phi2, Iphi1phi3;
		Iphi1phi2, Iphi2_2, Iphi2phi3;
		Iphi1phi3, Iphi2phi3, Iphi3_2
		];
	K_elem = [
		Idrv_phi1_2 , Idrv_phi1drv_phi2, Idrv_phi1drv_phi3;
		Idrv_phi1drv_phi2, Idrv_phi2_2, Idrv_phi2drv_phi3;
		Idrv_phi1drv_phi3, Idrv_phi2drv_phi3, Idrv_phi3_2
		];

	% global matrices construction
	M_global = sparse(2*N_elem+1,2*N_elem+1);
	K_global = sparse(2*N_elem+1,2*N_elem+1);
	for i=1:N_elem
		M_global(2*i-1:2*i+1,2*i-1:2*i+1) = M_global(2*i-1:2*i+1,2*i-1:2*i+1) + M_elem;
		K_global(2*i-1:2*i+1,2*i-1:2*i+1) = K_global(2*i-1:2*i+1,2*i-1:2*i+1) + K_elem;
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

	p = x(1:2*N_elem+1);
	R = x(2*N_elem+2);

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

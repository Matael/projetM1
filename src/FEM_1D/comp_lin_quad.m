%
% comp_lin_quad.m
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
R_vect_quad = zeros(N_elem_max,1);
R_vect_lin = zeros(N_elem_max,1);
for N_elem=1:N_elem_max

	% Discretization
	h = L/N_elem;

	%%%%%%%%%%%%%%%%%%%%%
	% BEGIN QUADRATIC
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
	M_elem_quad = [
		Iphi1_2 , Iphi1phi2, Iphi1phi3;
		Iphi1phi2, Iphi2_2, Iphi2phi3;
		Iphi1phi3, Iphi2phi3, Iphi3_2
		];
	K_elem_quad = [
		Idrv_phi1_2 , Idrv_phi1drv_phi2, Idrv_phi1drv_phi3;
		Idrv_phi1drv_phi2, Idrv_phi2_2, Idrv_phi2drv_phi3;
		Idrv_phi1drv_phi3, Idrv_phi2drv_phi3, Idrv_phi3_2
		];

	% global matrices construction
	M_global_quad = sparse(2*N_elem+1,2*N_elem+1);
	K_global_quad = sparse(2*N_elem+1,2*N_elem+1);
	for i=1:N_elem
		M_global_quad(2*i-1:2*i+1,2*i-1:2*i+1) = M_global_quad(2*i-1:2*i+1,2*i-1:2*i+1) + M_elem_quad;
		K_global_quad(2*i-1:2*i+1,2*i-1:2*i+1) = K_global_quad(2*i-1:2*i+1,2*i-1:2*i+1) + K_elem_quad;
	end

	% global matrix
	A_quad = sparse(2*N_elem+2,2*N_elem+2);
	A_quad(1:2*N_elem+1,1:2*N_elem+1) = k^2*M_global_quad - K_global_quad;

	A_quad(2*N_elem+2,1) = 1;
	A_quad(2*N_elem+2,2*N_elem+2) = -1;
	A_quad(1,2*N_elem+2) = -j*k;

	b_quad = zeros(2*N_elem+2,1);
	b_quad(1) = -j*k;
	b_quad(2*N_elem+2) = 1;
	%%%%%%%%%%%%%%%%%%%%%
	% END QUADRATIC
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%
	% BEGIN LINEAR
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	h = L/(2*N_elem);
	% elementary matrices
	M_elem_lin = h/6*[2 1; 1 2];
	K_elem_lin = 1/h*[1 -1; -1 1];

	% global matrices construction
	M_global_lin = sparse(2*N_elem+1,2*N_elem+1);
	K_global_lin = sparse(2*N_elem+1,2*N_elem+1);
	for i=1:2*N_elem
		M_global_lin(i:i+1,i:i+1) = M_global_lin(i:i+1,i:i+1) + M_elem_lin;
		K_global_lin(i:i+1,i:i+1) = K_global_lin(i:i+1,i:i+1) + K_elem_lin;
	end

	% global matrix
	A_lin = sparse(2*N_elem+2,2*N_elem+2);
	A_lin(1:2*N_elem+1,1:2*N_elem+1) = k^2*M_global_lin - K_global_lin;

	A_lin(2*N_elem+2,1) = 1;
	A_lin(2*N_elem+2,2*N_elem+2) = -1;
	A_lin(1,2*N_elem+2) = -j*k;

	b_lin = zeros(2*N_elem+2,1);
	b_lin(1) = -j*k;
	b_lin(2*N_elem+2) = 1;
	%%%%%%%%%%%%%%%%%%%%%
	% END LINEAR
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% resolution
	x_quad = A_quad\b_quad;
	R_vect_quad(N_elem) = x_quad(2*N_elem+2);
	x_lin = A_lin\b_lin;
	R_vect_lin(N_elem) = x_lin(2*N_elem+2);
end

% Analytical solution
Zi = Zc/(j*tan(k*L));
R_ana = (Zi-Zc)/(Zi+Zc);

% error
err_quad = zeros(N_elem_max, 1);
err_lin = zeros(N_elem_max, 1);
for i=1:N_elem_max
	err_quad(i) = norm(angle(R_vect_quad(i)) - angle(R_ana),2)/norm(angle(R_ana),2);
	err_lin(i) = norm(angle(R_vect_lin(i)) - angle(R_ana),2)/norm(angle(R_ana),2);
end


% figures
figure(1); % error (LogLog)
loglog(2*(1:N_elem_max) , err_quad, 'r', 'LineWidth', 2)
hold on;
loglog(2*(1:N_elem_max) , err_lin, 'b', 'LineWidth', 2)
xlabel('Degrees of freedom')
ylabel('Relative Error')
grid on;
set(gca, 'xminorgrid', 'off');
set(gca, 'yminorgrid', 'off');
legend('FEM Quad', 'FEM Lin')
print('-dpng', 'convergence.png');

figure(2); % R FEM & R ana (angle)
plot(2*(1:N_elem_max), angle(R_vect_quad)+(angle(R_vect_quad)<0)*2*pi, 'r+', 'LineWidth', 2);
hold on;
plot(2*(1:N_elem_max), angle(R_vect_lin)+(angle(R_vect_lin)<0)*2*pi, 'b+', 'LineWidth', 2);
plot([1 N_elem_max], [1 1]*angle(R_ana), 'k', 'LineWidth', 2);
xlabel('Degrees of freedom')
ylabel('Phase of R')
legend('FEM Quad', 'FEM Lin', 'Analytical')
xlim([0 125])
grid on;
print('-dpng', 'phase.png');

figure(3);
subplot(121); % error (LogLog)
loglog(2*(1:N_elem_max) , err_quad, 'r', 'LineWidth', 2)
hold on;
loglog(2*(1:N_elem_max) , err_lin, 'b', 'LineWidth', 2)
xlabel('Degrees of freedom')
ylabel('Relative Error')
grid on;
set(gca, 'xminorgrid', 'off');
set(gca, 'yminorgrid', 'off');
legend('FEM Quad', 'FEM Lin')

subplot(122); % R FEM & R ana (angle)
plot(2*(1:N_elem_max), angle(R_vect_quad)+(angle(R_vect_quad)<0)*2*pi, 'r+', 'LineWidth', 2);
hold on;
plot(2*(1:N_elem_max), angle(R_vect_lin)+(angle(R_vect_lin)<0)*2*pi, 'b+', 'LineWidth', 2);
plot([1 N_elem_max], [1 1]*angle(R_ana), 'k', 'LineWidth', 2);
xlabel('Degrees of freedom')
ylabel('Phase of R')
legend('FEM Quad', 'FEM Lin', 'Analytical')
xlim([0 125])
grid on;


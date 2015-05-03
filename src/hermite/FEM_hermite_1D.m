%
% FEM_1D_DGMlike_quad.m
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
addpath('../utils'); load_utils

% fluid params
rho = 1.2141;
c = 343;
Zc = rho*c;

% frequency
f = 500; % Hz
omega = 2*pi*f;
k = 2*pi*f/c;

% geometrical
L = 1;

% Analytical solution
Zi = Zc/(j*tan(k*L));
R_ana = (Zi-Zc)/(Zi+Zc);

N_elem_max=100;

err = zeros(N_elem_max, 1);
R_vect_ref = zeros(N_elem_max, 1);
R_vect_test = zeros(N_elem_max, 1);
err_R_ref = zeros(N_elem_max, 1);
err_R_test = zeros(N_elem_max, 1);

for N_elem=1:N_elem_max

	Nddl = 2*N_elem+1; % DDL quadratic FEM

	% Discretization
	h = L/N_elem;

	% form functions
	% %%%%%%%%%%%%%%%%%%%%%%%%% %
	%     QUADRATIC ELEMENTS    %
	% %%%%%%%%%%%%%%%%%%%%%%%%% %
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
	A = sparse(Nddl+1,Nddl+1);
	A(1:Nddl,1:Nddl) = k^2*M_global - K_global;

	A(Nddl+1,1) = 1;
	A(Nddl+1,Nddl+1) = -1;
	A(1,Nddl+1) = -j*k;

	b = zeros(Nddl+1,1);
	b(1) = -j*k;
	b(Nddl+1) = 1;

	% resolution
	x = A\b;

	p_ref = x(1:Nddl);
	R_vect_ref(N_elem) = x(Nddl+1);


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
	%        HERMITE POLYNOMIALS       %
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
	h = h/2; % we want to compare with the same number of DDL
	% form functions
	H0 = @(x) 2*(x/h)^3-3*(x/h)^2+1;
	H1 = @(x) ((x/h)^3-2*(x/h)^2+(x/h))*h;
	H2 = @(x) -2*(x/h)^3+3*(x/h)^2;
	H3 = @(x) ((x/h)^3-(x/h)^2)*h;


	drv_H0 = @(x) 6*x^2/h^3- 6*x/h^2;
	drv_H1 = @(x) (3*x^2/h^3-4^x/h^2+1)*h;
	drv_H2 = @(x) -6*x^2/h^3+6*x/h^2;
	drv_H3 = @(x) (3*x^2/h^3-2*x/h^2)*h;

	IH00 = intg_fun(mulfun(H0, H0), 0, h);
	IH01 = intg_fun(mulfun(H0, H1), 0, h);
	IH02 = intg_fun(mulfun(H0, H2), 0, h);
	IH03 = intg_fun(mulfun(H0, H3), 0, h);
	IH11 = intg_fun(mulfun(H1, H1), 0, h);
	IH12 = intg_fun(mulfun(H1, H2), 0, h);
	IH13 = intg_fun(mulfun(H1, H3), 0, h);
	IH22 = intg_fun(mulfun(H2, H2), 0, h);
	IH23 = intg_fun(mulfun(H2, H3), 0, h);
	IH33 = intg_fun(mulfun(H3, H3), 0, h);

	Idrv_H00 = intg_fun(mulfun(drv_H0, drv_H0), 0, h);
	Idrv_H01 = intg_fun(mulfun(drv_H0, drv_H1), 0, h);
	Idrv_H02 = intg_fun(mulfun(drv_H0, drv_H2), 0, h);
	Idrv_H03 = intg_fun(mulfun(drv_H0, drv_H3), 0, h);
	Idrv_H11 = intg_fun(mulfun(drv_H1, drv_H1), 0, h);
	Idrv_H12 = intg_fun(mulfun(drv_H1, drv_H2), 0, h);
	Idrv_H13 = intg_fun(mulfun(drv_H1, drv_H3), 0, h);
	Idrv_H22 = intg_fun(mulfun(drv_H2, drv_H2), 0, h);
	Idrv_H23 = intg_fun(mulfun(drv_H2, drv_H3), 0, h);
	Idrv_H33 = intg_fun(mulfun(drv_H3, drv_H3), 0, h);

	% elementary matrices
	M_elem = [
			IH00 IH01 IH02 IH03 ;
			IH01 IH11 IH12 IH13 ;
			IH02 IH12 IH22 IH23 ;
			IH03 IH13 IH23 IH33
		];
	K_elem = [
			Idrv_H00 Idrv_H01 Idrv_H02 Idrv_H03 ;
			Idrv_H01 Idrv_H11 Idrv_H12 Idrv_H13 ;
			Idrv_H02 Idrv_H12 Idrv_H22 Idrv_H23 ;
			Idrv_H03 Idrv_H13 Idrv_H23 Idrv_H33
		];

	% global matrices construction
	M_global = sparse(2*Nddl,2*Nddl);
	K_global = sparse(2*Nddl,2*Nddl);
	for i=1:N_elem*2
		M_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) = M_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) + M_elem;
		K_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) = K_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) + K_elem;
	end

	% global matrix
	A = k^2*M_global - K_global;

	b = zeros(2*Nddl,1);
	b(1) = - j*k;
	% A(1,1) = A(1,1) - 0.5*j*k*Zc*H0(0);
	% A(1,2) = A(1,2) - 0.5*j*k*H1(0)/(j*omega*rho);
	% A(1,3) = A(1,3) - 0.5*j*k*Zc*H2(0);
	% A(1,4) = A(1,4) - 0.5*j*k*H3(0)/(j*omega*rho);
	A(1,1) = A(1,1) - 0.5*j*k;
	A(1,2) = A(1,2) - 0.5/h;

	p_test = A\b;
	R_vect_test(N_elem) = Zc*0.5*(p_test(1)/Zc + p_test(2)/(j*omega*rho*h));

	% err(N_elem) = norm(p_ref-p_test, 2)/norm(p_ref,2);
	err_R_ref(N_elem) = norm(angle(R_vect_ref(N_elem))-angle(R_ana),2)/norm(angle(R_ana),2);
	err_R_test(N_elem) = norm(angle(R_vect_test(N_elem))-angle(R_ana),2)/norm(angle(R_ana),2);
end

% figure;
% plot(real(p_ref), '+r', 'LineWidth', 2);
% hold on;
% plot(real(p_test), 'ob', 'LineWidth', 2);

% figure;
% loglog(1:N_elem_max, err, '+r', 'LineWidth', 2)
% xlabel('Nb DOF')
% grid on;

figure;
loglog(1:N_elem_max, err_R_ref, 'b', 'LineWidth', 2)
hold on;
loglog(1:N_elem_max, err_R_test, 'r', 'LineWidth', 2)
legend('Quad', 'Hermite')
grid on;
xlabel('Nb. Elements')
ylabel('Error')
% print('-dpng', 'comparison_DGMlike_FEM.png')

figure;
plot(real(p_test(2*(1:Nddl)-1)), 'LineWidth', 2)
hold on;
plot(real(p_ref), 'r')

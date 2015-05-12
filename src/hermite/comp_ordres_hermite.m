%
% FEM_1D_hermite.m
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

% frequency settings
f = 500;% Hz

% fluid settings
rho = 1.241;
c = 343;
Zc = rho*c;

omega = 2*pi*f;
k = omega/c;

% Geometrical settings
L = 1;

% Analytical solution
Zi = Zc/(j*tan(k*L));
R_ana = (Zi-Zc)/(Zi+Zc);

N_elem_max=100;

R_vect_ref = zeros(1,N_elem_max);
R_vect_test = zeros(1,N_elem_max);

for N_elem=1:N_elem_max

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%                       H-Splines Interpolation                      %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	N_elem_herm = 2*N_elem-1; % We want the same number of ddl
	Nddl_herm = N_elem_herm+1;

	h = L/N_elem_herm;
	% form functions
	H0 = @(x) 2*(x/h)^3-3*(x/h)^2+1;
	H1 = @(x) ((x/h)^3-2*(x/h)^2+(x/h))*h;
	H2 = @(x) -2*(x/h)^3+3*(x/h)^2;
	H3 = @(x) ((x/h)^3-(x/h)^2)*h;


	drv_H0 = @(x) 6*x^2/h^3- 6*x/h^2;
	drv_H1 = @(x) (3*x^2/h^3-4*x/h^2+1/h)*h;
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


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%      USING CHARACTERISTICS FORMULATION      %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% global matrices construction
	M_global = sparse(2*Nddl_herm,2*Nddl_herm);
	K_global = sparse(2*Nddl_herm,2*Nddl_herm);
	for i=1:N_elem_herm
		% debug stuff
		% msg = sprintf('iter[%d] %d -> %d', i, 2*i-1, 2*(i+1));
		% disp(msg)
		%%%%%%%%%
		M_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) = M_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) + M_elem;
		K_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) = K_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) + K_elem;
	end

	% global matrix
	A = k^2*M_global - K_global;

	b = zeros(2*Nddl_herm,1);
	b(1) = -j*k;
	A(1,1) = A(1,1) - 0.5*j*k;
	A(1,2) = A(1,2) - 0.5;

	p_test = A\b;
	R_vect_test(N_elem) = 0.5*(p_test(1)+ p_test(2)/(j*k));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%    USING CLASSICAL FEM W/ R FORMULATION     %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% global matrices construction
	M_global = sparse(2*Nddl_herm,2*Nddl_herm);
	K_global = sparse(2*Nddl_herm,2*Nddl_herm);
	for i=1:N_elem_herm
		M_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) = M_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) + M_elem;
		K_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) = K_global(2*i-1:2*(i+1),2*i-1:2*(i+1)) + K_elem;
	end

	% global matrix
	A = sparse(2*Nddl_herm+1,2*Nddl_herm+1);
	A(1:2*Nddl_herm,1:2*Nddl_herm) = k^2*M_global - K_global;

	b = zeros(2*Nddl_herm+1,1);
	b(1) = -j*k;
	b(2*Nddl_herm+1) = 1;
	A(1,2*Nddl_herm+1) = -j*k;
	A(2*Nddl_herm+1,1) = 1;
	A(2*Nddl_herm+1,2*Nddl_herm+1) = -1;

	p_ref = A\b;
	R_vect_ref(N_elem) = p_ref(2*Nddl_herm+1);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%                           Errors Computation                          %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	err_R_ref(N_elem) = norm(angle(R_vect_ref(N_elem))-angle(R_ana),2)/norm(angle(R_ana),2);
	err_R_test(N_elem) = norm(angle(R_vect_test(N_elem))-angle(R_ana),2)/norm(angle(R_ana),2);
end


absisses = 2*(1:N_elem_max)+1;

figure;
loglog(absisses,err_R_ref, 'b', 'LineWidth', 2)
hold on;
loglog(absisses,err_R_test, 'r', 'LineWidth', 2)
grid on;
xlabel('Degres of Freedom')
ylabel('log(relative error)')
legend('FEM+R', 'FEM+Characteristics');
print('-dpng', 'comparaison_ordres_hermite.png')

% figure;
% plot(absisses,angle(R_vect_ref), 'b', 'LineWidth', 2)
% hold on;
% plot(absisses,angle(R_vect_test), 'r', 'LineWidth', 2)
% grid on;
% xlabel('Degres of Freedom')
% ylabel('Angle(R)')

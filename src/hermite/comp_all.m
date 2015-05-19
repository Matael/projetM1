%
% comp_all.m
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

comp_ordres_hermite

herm_pR = err_R_ref;
herm_car = err_R_test;

figure(42);
elems = 1:N_elem_max;
hold on;
loglog(elems, herm_car, 'r', 'LineWidth', 2)
loglog(elems, herm_pR, '--b', 'LineWidth', 2)

addpath('../FEM_DGM_form/')
FEM_1D_DGMlike_quad;
quad_pR =  err_R_ref;
quad_car = err_R_test;

figure(42)
elems = 1:N_elem_max;
loglog(elems, quad_pR, '--k', 'LineWidth', 2)
loglog(elems, quad_car, 'k', 'LineWidth', 2)
xlabel('Degrees of freedom')
ylabel('Relative Error')
grid on;
print('-dpng', 'herm_comp.png')






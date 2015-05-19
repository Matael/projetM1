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

% clear all;
% close all;

FEM_1D_DGMlike_quad;

quad_errref = err_R_ref;
quad_errtest = err_R_test;

figure;
elems = 1:N_elem_max;
hold on;
loglog(elems, quad_errref, 'k', 'LineWidth', 2)
loglog(elems, quad_errtest, 'r', 'LineWidth', 2)
xlabel('Degrees of freedom')
ylabel('Relative Error')
grid on;
print('-dpng', 'FEM_DGMlike_quad.png')



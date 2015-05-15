%
% comp_FEMherm_dgm.m
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

addpath('../hermite');
DGM_1D_PW;
FEM_1D_hermite;

p_test = real(p_test((1:Nddl_herm)*2-1))'; % FEM
p = real(p); % DGM

figure;
plot(x_v,p_test,'b');
hold on;
plot(x_v,p,'r');
legend('Hermite FEM', 'DGM')
print('-dpng', 'comp_hermiteFEM_dgm.png')


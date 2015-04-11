%
% intg_fun.m
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

function [I] = intg_fun(fun,a,b)
	n_pg = 4;

	points_gauss(1)= 0.339981043584856;
	points_gauss(2)=-0.339981043584856;
	points_gauss(3)= 0.861136311594053;
	points_gauss(4)=-0.861136311594053;

	weight_gauss(1)=0.652145154862546;
	weight_gauss(2)=0.652145154862546;
	weight_gauss(3)=0.347854845137454;
	weight_gauss(4)=0.347854845137454;

	I = 0;
	h = b-a;
	for i_pg=1:n_pg
		xi=points_gauss(i_pg);
		I = I + fun((b-a)/2*xi+(a+b)/2)*weight_gauss(i_pg);
	end

	I = I*h/2;
end




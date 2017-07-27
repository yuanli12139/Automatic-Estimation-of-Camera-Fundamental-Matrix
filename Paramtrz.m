% Parameterization of homogenous vectors %

function [ v ] = Paramtrz( v_bar )

v_bar_T=v_bar';
v_bar_T=v_bar_T/norm(v_bar_T);
a=v_bar_T(1);
b=v_bar_T(2:end)';

v=2/sinc(acos(a)/pi)*b;

vn=norm(v,'fro');
if vn>pi
    v=(1-2*pi/vn*ceil((vn-pi)/(2*pi)))*v;
end

end


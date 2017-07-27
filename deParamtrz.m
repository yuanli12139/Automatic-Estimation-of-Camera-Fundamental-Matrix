% Deparameterization of homogenous vectors %

function [ v_bar ] = deParamtrz( v )

vn=norm(v,'fro');
a=cos(vn/2);
b=sinc(vn/2/pi)/2*v;

v_bar=[a,b']';

end


function [ X_pi ] = getXpi( xh,Pi )

a=Pi(1);b=Pi(2);c=Pi(3);d=Pi(4);

x=xh(1);y=xh(2);w=xh(3);

X_pi=[d*x;d*y;d*w;-(a*x+b*y+c*w)];

end


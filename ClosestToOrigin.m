function [ x ] = ClosestToOrigin( l )

a=l(1);
b=l(2);
c=l(3);

x=[-a*c;-b*c;a^2+b^2];

end


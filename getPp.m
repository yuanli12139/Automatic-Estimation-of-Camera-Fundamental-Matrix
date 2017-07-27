function [ Pp ] = getPp( F )

[U,D,V]=svd(F);
diagD=diag(D);
s=diagD(1);
t=diagD(2);

% permutation matrix
W=[0,1,0;
   -1,0,0;
   0,0,0];

Z=[0,-1,0;
   1,0,0;
   0,0,1];

Dp=diag([s,t,(s+t)/2]');

S=U*W*U';
M=U*Z*Dp*V';

a1=S(3,2);a2=S(1,3);a3=S(2,1);
ep=[a1;a2;a3];

Pp=[M,ep];

end


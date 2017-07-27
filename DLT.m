function [ F_DLT ] = DLT( x1,x12,numPts ) %Input: inhomoPts (inliers)

% Data Normalization
[xp,T]=Norml(x12);
[x,You]=Norml(x1);

A=zeros(numPts,9); % 8+ points
for pt=1:numPts
    A(pt,:)=kron(xp(:,pt)',x(:,pt)');
end

[U,D,V]=svd(A);
V_T=V';
f=V_T(end,:);

FMatrix=reshape(f,3,3)';

[U_F,D_F,V_F]=svd(FMatrix);
Dp=D_F;
Dp(3,3)=0;
FMatrix=U_F*Dp*V_F';

F_DLT_unscaled=deNorml(FMatrix,T,You);    % denormalization
F_DLT=F_DLT_unscaled/norm(F_DLT_unscaled,'fro');

format longg;
F_DLT

end
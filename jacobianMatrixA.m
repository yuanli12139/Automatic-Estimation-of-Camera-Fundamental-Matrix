% ------------- Jacobian Matrix ------------ %
% X: 3D homogenous vectors
% v: deparameterized (n-1)-vector
% P: 3x4 projection matrix

function [ J ] = jacobianMatrixA( Xh, P )

X=deParamtrz(Xh);
P_T=P';v=Paramtrz(P_T(:));

x_ho=P*X;
x_inho=[x_ho(1,:)./x_ho(3,:);x_ho(2,:)./x_ho(3,:)];

n=size(x_inho,2);
J=[];
for i=1:n

% 2 by 12 matrix:
    w=P(3,1)*X(1,i)+P(3,2)*X(2,i)+P(3,3)*X(3,i)+P(3,4)*X(4,i);
    pxh_pvb=[1/w*X(:,i)',zeros(4,1)',-x_inho(1,i)/w*X(:,i)';zeros(4,1)',1/w*X(:,i)',-x_inho(2,i)/w*X(:,i)'];

% 12 by 11 matrix:   
    vn=norm(v,'fro');
    da_pv=-sinc(vn/2/pi)/4*v';
    
    if vn==0   
       pb_pv=1/2*eye(length(v)); 
    else
       pb_pv=sinc(vn/2/pi)/2*eye(length(v))+1/(4*vn)*(cos(vn/2)/(vn/2)-sin(vn/2)/(vn/2)^2)*(v*v');    
    end
    
    pvb_pv=[da_pv;pb_pv];
    
% Jacobian matrix (2n by 11):    
    J=[J;pxh_pvb*pvb_pv];
end

end


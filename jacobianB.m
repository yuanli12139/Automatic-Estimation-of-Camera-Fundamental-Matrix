function [ J ] = jacobianB( Xh,P )

% 2 by 4 matrix:
    Xb=deParamtrz(Xh);
    %xh=homo2inhomo(P*Xb);

    syms X Y Z T p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12
    Partial1=jacobian([(X*p1+Y*p2+Z*p3+T*p4)/(X*p9+Y*p10+Z*p11+T*p12),(X*p5+Y*p6+Z*p7+T*p8)/(X*p9+Y*p10+Z*p11+T*p12)],[X,Y,Z,T]);
    pXh_pXb=subs(Partial1,[X Y Z T p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12],[Xb(1),Xb(2),Xb(3),Xb(4),P(1,1),P(1,2),P(1,3),P(1,4),P(2,1),P(2,2),P(2,3),P(2,4),P(3,1),P(3,2),P(3,3),P(3,4)]);
    pXh_pXb=double(pXh_pXb);
    
% 4 by 3 matrix:   
    Xhn=norm(Xh,'fro');
    da_pXh=-sinc(Xhn/2/pi)/4*Xh';
    
    if Xhn==0   
       pb_pXh=1/2*eye(length(Xh));
    else
       pb_pXh=sinc(Xhn/2/pi)/2*eye(length(Xh))+1/(4*Xhn)*(cos(Xhn/2)/(Xhn/2)-sin(Xhn/2)/(Xhn/2)^2)*(Xh*Xh');    
    end
    
    pXb_pXh=[da_pXh;pb_pXh];
    
% Jacobian matrix (2 x 3): 
    J=pXh_pXb*pXb_pXh;
    
end


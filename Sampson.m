function [ SPerror,dlt ] = Sampson( xMatch1,yMatch1,xMatch12,yMatch12,F )

    n=numel(xMatch12); 
    SPerror=zeros(1,n);  
     
    x=inhomo2homo([xMatch1';yMatch1']);
    xp=inhomo2homo([xMatch12';yMatch12']);
    
    for i=1:n
        ep=xp(:,i)'*F*x(:,i);
        
        xi=xMatch1(i);yi=yMatch1(i);
        xpi=xMatch12(i);ypi=yMatch12(i);
        J=[xpi*F(1,1)+ypi*F(2,1)+F(3,1),xpi*F(1,2)+ypi*F(2,2)+F(3,2),xi*F(1,1)+yi*F(1,2)+F(1,3),xi*F(2,1)+yi*F(2,2)+F(2,3)];
        
        lbd=-ep/(J*J');
        dlt(:,i)=J'*lbd;

        SPerror(i)=dlt(:,i)'*dlt(:,i);
    end

end


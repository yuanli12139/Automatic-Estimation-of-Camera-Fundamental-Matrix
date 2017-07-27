function [ X_pi ] = Triangulation( xIn1,yIn1,xIn12,yIn12,F,Pp )

xs=inhomo2homo([xIn1';yIn1']);
xps=inhomo2homo([xIn12';yIn12']);

n=size(xs,2);

X_pi=zeros(4,n);

for pt=1:n
    x=xs(1,pt);y=xs(2,pt);w=xs(3,pt);
    xp=xps(1,pt);yp=xps(2,pt);wp=xps(3,pt);
    
    T=[w,0,-x;0,w,-y;0,0,w];
    Tp=[wp,0,-xp;0,wp,-yp;0,0,wp];
    
    Fs=inv(Tp')*F*inv(T);
    
    e=null(Fs);e1=e(1);e2=e(2);
    ep=null(Fs');ep1=ep(1);ep2=ep(2);

    e=sqrt(1/(e1^2+e2^2))*e;e1=e(1);e2=e(2);e3=e(3);
    ep=sqrt(1/(ep1^2+ep2^2))*ep;ep1=ep(1);ep2=ep(2);ep3=ep(3);
    
    R=[e1,e2,0;-e2,e1,0;0,0,1];
    Rp=[ep1,ep2,0;-ep2,ep1,0;0,0,1];
    
    Fs=Rp*Fs*R';
    
    a=Fs(2,2);b=Fs(2,3);c=Fs(3,2);d=Fs(3,3);f=e3;fp=ep3;
    
    syms as bs c_s ds fs fps ts
    g_t=ts*((as*ts+bs)^2+fps^2*(c_s*ts+ds)^2)^2-(as*ds-bs*c_s)*(1+fs^2*ts^2)^2*(as*ts+bs)*(c_s*ts+ds);
    
    polyn=collect(g_t,ts);
    
    Eq=subs(polyn,[as,bs,c_s,ds,fs,fps],[a,b,c,d,f,fp]);
    [cefs,~]=coeffs(Eq,ts);
    t=real(roots(double(cefs)));
    
    s_t=t.^2./(1+f^2.*t.^2)+(c*t+d).^2./((a*t+b).^2+fp^2*(c*t+d).^2);
    t=t(find(s_t==min(s_t)));
    
    l=[t*f,1,-t]';
    lp=[-fp*(c*t+d),a*t+b,c*t+d]';
    
    xh=ClosestToOrigin(l);
    xph=ClosestToOrigin(lp);
    
    xh=inv(T)*R'*xh;
    xph=inv(Tp)*Rp'*xph;
    
    l_perp=getLperp(xh,xph,F);
    
    Pi=Pp'*l_perp;
    
    X_pi(:,pt)=getXpi(xh,Pi);
    
end

end


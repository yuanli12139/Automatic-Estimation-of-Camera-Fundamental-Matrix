function [ F_final,numSolutions ] = F7( x1,x2,x3,x4,x5,x6,x7 , x12,x22,x32,x42,x52,x62,x72 )
%function [ alf_real,numSolutions ] = F7( x1,x2,x3,x4,x5,x6,x7 , x12,x22,x32,x42,x52,x62,x72 )
x=[x1,x2,x3,x4,x5,x6,x7];
xp=[x12,x22,x32,x42,x52,x62,x72];

A=zeros(7,9);
for pt=1:7
    A(pt,:)=kron(xp(:,pt)',x(:,pt)');
end

f=null(A);
a=f(:,1);b=f(:,2);

syms a1 a2 a3 a4 a5 a6 a7 a8 a9
F1 = [ a1, a2, a3; a4, a5, a6; a7, a8, a9 ];
syms b1 b2 b3 b4 b5 b6 b7 b8 b9
F2 = [ b1, b2, b3; b4, b5, b6; b7, b8, b9 ];
syms alf
F = alf * F1 + F2;
g = det( F );
polyn=collect( g, alf );

Eq=subs(polyn,[a1,a2,a3,a4,a5,a6,a7,a8,a9,b1,b2,b3,b4,b5,b6,b7,b8,b9],[a',b']);
[c,~]=coeffs(Eq);

r=roots(double(c));
tol=1e-100;
alf_real = real(r(tol>abs(imag(r))));
%alf_real=alf(find(imag(alf)==0));

numSolutions=numel(alf_real);

F1=reshape(a,3,3)';
F2=reshape(b,3,3)';

F_final=cell(1,numSolutions);
for i=1:numel(alf_real)
    F_final{1,i}=(alf_real(i)*F1+F2)/norm((alf_real(i)*F1+F2),'fro');
end

end


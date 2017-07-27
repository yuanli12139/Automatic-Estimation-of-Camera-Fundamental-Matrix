function [ minorEigenImg ] = minorEigen( Ix,Iy,winSize )

minorEigenImg=zeros(size(Ix));
p=floor(winSize/2);
Ix_pad=padarray(Ix,[p,p]);
Iy_pad=padarray(Iy,[p,p]);
for row=p+1:size(Ix,1)+p
    for col=p+1:size(Ix,2)+p
        E_Ix2=sum(sum(Ix_pad(row-p:row+p,col-p:col+p).^2));
        E_Iy2=sum(sum(Iy_pad(row-p:row+p,col-p:col+p).^2));
        E_IxIy=sum(sum(Ix_pad(row-p:row+p,col-p:col+p).*Iy_pad(row-p:row+p,col-p:col+p)));
        N=[E_Ix2,E_IxIy;
           E_IxIy,E_Iy2]/winSize^2;
        
        minorEigenImg(row-p,col-p)=min(eig(N));
    end
end

end


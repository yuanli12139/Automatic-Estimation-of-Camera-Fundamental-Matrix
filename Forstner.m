function [ xCorner,yCorner ] = Forstner( img_nms,Ix,Iy,winSize )

xCorner=[];
yCorner=[];

p=floor(winSize/2);

img_nmsPd=padarray(img_nms,[p,p]);
IxPd=padarray(Ix,[p,p]);
IyPd=padarray(Iy,[p,p]);

[row,col]=find(img_nmsPd~=0);

l=size(Ix,2);
w=size(Ix,1);

[xCor,yCor]=meshgrid(1:l,1:w);
xCorPd=padarray(xCor,[p,p]);
yCorPd=padarray(yCor,[p,p]);

for fp=1:size(row,1)
    
    Ix_w=IxPd(row(fp)-p:row(fp)+p,col(fp)-p:col(fp)+p);
    Iy_w=IyPd(row(fp)-p:row(fp)+p,col(fp)-p:col(fp)+p);
    
    xCor_w=xCorPd(row(fp)-p:row(fp)+p,col(fp)-p:col(fp)+p);
    yCor_w=yCorPd(row(fp)-p:row(fp)+p,col(fp)-p:col(fp)+p);
    
    b=[sum(sum(xCor_w.*Ix_w.^2+yCor_w.*Ix_w.*Iy_w));
       sum(sum(xCor_w.*Ix_w.*Iy_w+yCor_w.*Iy_w.^2))]/winSize^2;
    
    A=[sum(sum(Ix_w.^2)),sum(sum(Ix_w.*Iy_w));
       sum(sum(Ix_w.*Iy_w)),sum(sum(Iy_w.^2))]/winSize^2;
    
    Corner=A\b;
    xCorner=[xCorner;Corner(1)];
    yCorner=[yCorner;Corner(2)];
    
end

end


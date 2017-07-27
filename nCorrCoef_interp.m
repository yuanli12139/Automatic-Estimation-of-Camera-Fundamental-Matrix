function [ xMatch1,yMatch1,xMatch12,yMatch12 ] = nCorrCoef_interp( img1,img2,x1,y1,x2,y2,Thres,winSize )

rows1=floor(y1);rows2=floor(y2);
cols1=floor(x1);cols2=floor(x2);
p=floor(winSize/2);

[Xq,Yq]=meshgrid(1:1:1+winSize);

ptsCount=0;

corrs12=zeros(numel(rows2),numel(rows1));
for row1=1:numel(rows1)
    for row2=1:numel(rows2)
        win1=img1(rows1(row1)-p:rows1(row1)+p+1,cols1(row1)-p:cols1(row1)+p+1);
        X1=meshgrid((y1(row1)-rows1(row1)+1:1:y1(row1)-rows1(row1)+winSize));
        Y1=meshgrid((y1(row1)-rows1(row1)+1:1:y1(row1)-rows1(row1)+winSize))';
        win1_ip=interp2(Xq,Yq,win1,X1,Y1);
        
        win2=img2(rows2(row2)-p:rows2(row2)+p+1,cols2(row2)-p:cols2(row2)+p+1);
        X2=meshgrid((y2(row2)-rows2(row2)+1:1:y2(row2)-rows2(row2)+winSize));
        Y2=meshgrid((y2(row2)-rows2(row2)+1:1:y2(row2)-rows2(row2)+winSize))';
        win2_ip=interp2(Xq,Yq,win2,X2,Y2);
        
        corrs12(row2,row1)=corr2(win1_ip,win2_ip);
    end
    ptsCount=ptsCount+1
end

%corrs21=corrs12';

%{
wins1=zeros(winSize^2,numel(rows1));
wins2=zeros(winSize^2,numel(rows2));
%s=struct
for i=1:numel(rows1)
    win1=img1(rows1(i)-p:rows1(i)+p,cols1(i)-p:cols1(i)+p);
    wins1(:,i)=win1(:);
end

for i=1:numel(rows2)
    win2=img2(rows2(i)-p:rows2(i)+p,cols2(i)-p:cols2(i)+p);
    wins2(:,i)=win2(:);
end

corrs12=zeros(numel(rows2),numel(rows1));
for i=1:numel(rows1)
    corrs12(:,i)=arrayfun(@(index) sum(abs(normxcorr2(wins1(:,i),wins2(:,index)))),1:size(wins2,2));
end
toc;

%corrs12t=arrayfun(@(index) sum(normxcorr2(wins1,wins2(:,index))),1:size(wins2,2));

corrs21=zeros(numel(rows1),numel(rows2));
for i=1:numel(rows2)
    corrs21(:,i)=arrayfun(@(index) sum(abs(normxcorr2(wins2(:,i),wins1(:,index)))),1:size(wins1,2));
end
toc;
%}
maxCorrs12=max(corrs12);
idxMaxCorrs12=arrayfun(@(index) find(corrs12(:,index)==maxCorrs12(index),1),1:size(corrs12,2));
f1=find(maxCorrs12>=Thres);
fMatch12=idxMaxCorrs12(f1);
girls=zeros(3,numel(x2));
girls(1,:)=1:numel(x2);
for n=1:numel(fMatch12)
    if maxCorrs12(f1(n))>girls(3,fMatch12(n));
        girls(2,fMatch12(n))=f1(n);
        girls(3,fMatch12(n))=maxCorrs12(f1(n));
    end
end

girls2=girls(:,girls(3,:)>0);
couple=girls2(1:2,:);

xMatch1=x1(couple(2,:));yMatch1=y1(couple(2,:));
xMatch12=x2(couple(1,:));yMatch12=y2(couple(1,:));

%{
xMatch1=x1(f1);yMatch1=y1(f1);
xMatch12=x2(fMatch12);yMatch12=y2(fMatch12);

maxCorrs21=max(corrs21);
idxMaxCorrs21=arrayfun(@(index) find(corrs21(:,index)==maxCorrs21(index),1),1:size(corrs21,2));
f2=find(maxCorrs21>=Thres);
fMatch21=idxMaxCorrs21(f2);
xMatch2=x2(f2);yMatch2=y2(f2);
xMatch21=x1(fMatch21);yMatch21=y1(fMatch21);
%}
end
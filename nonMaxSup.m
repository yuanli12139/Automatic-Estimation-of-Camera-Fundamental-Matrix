function [ img_nms,count ] = nonMaxSup( img,Thres,winSize )

p=floor(winSize/2);
img_nms=zeros(size(img));
img_pad=padarray(img,[p,p]);
for row=p+1:size(img,1)+p
    for col=p+1:size(img,2)+p
        win=img_pad(row-p:row+p,col-p:col+p);
        [r,c]=find(win==max(max(win)));
        img_pad(row-p:row+p,col-p:col+p)=0;
        if max(max(win))>=Thres  
            for pt=1:numel(r)
                img_pad(r(pt)+row-p-1,c(pt)+col-p-1)=max(max(win));
            end
        end
    end
end

img_nms=img_pad(p+1:size(img_nms,1)+p,p+1:size(img_nms,2)+p);
count=sum(sum(img_nms~=0));

end
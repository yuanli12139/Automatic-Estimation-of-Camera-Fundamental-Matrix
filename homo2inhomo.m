function [ inhomo ] = homo2inhomo( homo )

row=size(homo,1);
col=size(homo,2);
inhomo=zeros(row-1,col);
for l=1:row-1

    inhomo(l,:)=homo(l,:)./homo(end,:);
    
end


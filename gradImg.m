function [ Ix,Iy ] = gradImg ( img )

k=1/12*[-1,8,0,-8,1]';  %Forstner corner point operator

Ix=conv2(img,k','same');
Iy=conv2(img',k','same');Iy=Iy';

end

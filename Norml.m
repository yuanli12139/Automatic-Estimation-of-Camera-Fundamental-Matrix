%----------Data Normalization----------%

function [ x_normalized, N ] = Norml( inhomoPts )

if size(inhomoPts,2) == 2
    
    mu2=mean(inhomoPts,1);
    mux2=mu2(:,1);
    muy2=mu2(:,2);

    sig2=sqrt(sum(var(inhomoPts,1)));
    s2=sqrt(2)/sig2;

    N=[s2,0,-mux2*s2;0,s2,-muy2*s2;0,0,1];
    
    homoPts=[inhomoPts,ones(size(inhomoPts,1),1)];
    x_normalized=N*homoPts';
    
end

if size(inhomoPts,2) == 3
    mu3=mean(inhomoPts,1);
    mux3=mu3(:,1);
    muy3=mu3(:,2);
    muz3=mu3(:,3);
    sig3=sqrt(sum(var(inhomoPts,1)));
    s3=sqrt(3)/sig3;

    N=[s3,0,0,-mux3*s3;0,s3,0,-muy3*s3;0,0,s3,-muz3*s3;0,0,0,1];

    homoPts=[inhomoPts,ones(size(inhomoPts,1),1)];
    x_normalized=N*homoPts';
end
end
function [ pk3Pts ] = pickOutliers( xCorner1_rb,yCorner1_rb,xIn1,yIn1 )
% pick 3 outliers
checkInliers=1;
while checkInliers
    pk3Pts=randperm(numel(xCorner1_rb),3);
    x_test=[xCorner1_rb(pk3Pts),yCorner1_rb(pk3Pts)]';
    for i=1:3
        checkInliers=find(x_test(1,i)==xIn1);
        if ~isempty(checkInliers)
            checkInliers=x_test(2,i)==yIn1(checkInliers);
        end
    end
end

end
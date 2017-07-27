function [ xIn1,yIn1,xIn12,yIn12,numInliers,Trials,dltIn ] = MSAC( xMatch1,yMatch1,xMatch12,yMatch12 )

n=numel(xMatch12);

S=7;  %Sample size
p=0.99;  %Probability that at least one of the random samples does not contain any outliers 
alf=0.95;  %Probability that the data point is an inlier
sigma2=1;  %Variance of the measurement error
m=1;  %Codimension of the model
Tolerence=chi2inv(alf,m);  %The square distance threshold
consensusMinCost=Inf;  %Initialize the consensusMinCost
maxTrials=Inf;  %Initialize the maxTrials
Threshold=0;

Trials=0;
while (Trials < maxTrials) && (consensusMinCost > Threshold) 
    Trials=Trials+1;
    if mod(Trials,20)==0 
        Trials 
        maxTrials
    end
    
    %Select a random sample (7-point algorithm)
    pickPts=randperm(n,S);
    x1=inhomo2homo([xMatch1(pickPts(1));yMatch1(pickPts(1))]);
    x2=inhomo2homo([xMatch1(pickPts(2));yMatch1(pickPts(2))]);
    x3=inhomo2homo([xMatch1(pickPts(3));yMatch1(pickPts(3))]);
    x4=inhomo2homo([xMatch1(pickPts(4));yMatch1(pickPts(4))]);
    x5=inhomo2homo([xMatch1(pickPts(5));yMatch1(pickPts(5))]);
    x6=inhomo2homo([xMatch1(pickPts(6));yMatch1(pickPts(6))]);
    x7=inhomo2homo([xMatch1(pickPts(7));yMatch1(pickPts(7))]);
    
    x12=inhomo2homo([xMatch12(pickPts(1));yMatch12(pickPts(1))]);
    x22=inhomo2homo([xMatch12(pickPts(2));yMatch12(pickPts(2))]);
    x32=inhomo2homo([xMatch12(pickPts(3));yMatch12(pickPts(3))]);
    x42=inhomo2homo([xMatch12(pickPts(4));yMatch12(pickPts(4))]);
    x52=inhomo2homo([xMatch12(pickPts(5));yMatch12(pickPts(5))]);
    x62=inhomo2homo([xMatch12(pickPts(6));yMatch12(pickPts(6))]);
    x72=inhomo2homo([xMatch12(pickPts(7));yMatch12(pickPts(7))]);
  
    %Calculate the model: Fundamental matrix (7-point algorithm)
    [F,numSolutions]=F7(x1,x2,x3,x4,x5,x6,x7 , x12,x22,x32,x42,x52,x62,x72);
    
    %Calculate error each datapoint (Sampson error)
    if numSolutions==3 
        SPerror=zeros(3,n);
        dlt=cell(1,3);
        [ SPerror(1,:),dlt{1} ] = Sampson( xMatch1,yMatch1,xMatch12,yMatch12,F{1} );
        [ SPerror(2,:),dlt{2} ] = Sampson( xMatch1,yMatch1,xMatch12,yMatch12,F{2} );
        [ SPerror(3,:),dlt{3} ] = Sampson( xMatch1,yMatch1,xMatch12,yMatch12,F{3} );
        %Calculate the cost (MSAC) 
        cost0=0;
        cst=zeros(1,3);
        for s=1:3
            for pt=1:n
                if SPerror(s,pt)<=Tolerence
                    cost0=cost0+SPerror(s,pt);
                else
                    cost0=cost0+Tolerence;
                end
                cst(pt)=cost0;
            end
        end
        minIdx=find(cst==min(cst),1);
        dlt=dlt{minIdx};
        cst=min(cst);
    else
        SPerror=zeros(1,n);
        dlt=zeros(4,n);
        [ SPerror,dlt ] = Sampson( xMatch1,yMatch1,xMatch12,yMatch12,F{1} );
        %Calculate the cost (MSAC) 
        cost0=0;
        for pt=1:n
            if SPerror(pt)<=Tolerence
                cost0=cost0+SPerror(pt);
            else
                cost0=cost0+Tolerence;
            end            
        end
        minIdx=1;
        cst=cost0;
    end
   
    if cst < consensusMinCost
        consensusMinCost=cst;
        consensusMinModel=F{minIdx};
        Inliers=find(SPerror(minIdx,:)<=Tolerence);
        numInliers=numel(find(SPerror(minIdx,:)<=Tolerence));

        xIn1=xMatch1(Inliers);yIn1=yMatch1(Inliers);
        xIn12=xMatch12(Inliers);yIn12=yMatch12(Inliers);
        dltIn=dlt(:,Inliers);

        w=numInliers/n;
        maxTrials=log10(1-p)/log10(1-w^S);
    end
end

end


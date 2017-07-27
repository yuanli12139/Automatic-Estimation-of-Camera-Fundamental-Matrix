% ---- Automatic estimation of the fundamental matrix ---- %
clear all;close all;
tic;

format longg;
%% (a) Feature detection %%
img1=imread('IMG_5030.JPG');img2=imread('IMG_5031.JPG');
img1=double(rgb2gray(img1));img2=double(rgb2gray(img2));

figure(1);
subplot(121);imshow(uint8(img1));
title('Input Image #1','FontSize',16);
subplot(122);imshow(uint8(img2));
title('Input Image #2','FontSize',16);

%[I1x,I1y]=gradImg(img1);
%[I2x,I2y]=gradImg(img2);
load('gradImg');
winSize_me=9;
%img1_me=minorEigen(I1x,I1y,winSize_me);
%img2_me=minorEigen(I2x,I2y,winSize_me);
load('minorEigen');
Thres_nms=10.7;
winSize_nms=9;
%[img1_nms,count1]=nonMaxSup(img1_me,Thres_nms,winSize_nms);
%[img2_nms,count2]=nonMaxSup(img2_me,Thres_nms,winSize_nms);
%save('nonMaxSup_new.mat','img1_nms','count1','img2_nms','count2');
load('nonMaxSup_new');
%[xCorner1,yCorner1]=Forstner(img1_nms,I1x,I1y,winSize_nms);
%[xCorner2,yCorner2]=Forstner(img2_nms,I2x,I2y,winSize_nms);
%save('Forstner_new.mat','xCorner1','yCorner1','xCorner2','yCorner2');
load('Forstner_new');

figure(2);
subplot(121);imshow(uint8(img1));
title('Image #1','FontSize',16);hold on;
for n=1:size(xCorner1,1)   
    rectangle('Position',[xCorner1(n)-4,yCorner1(n)-4,9,9],'LineWidth',1.1,'EdgeColor','y');
end
hold off;
subplot(122);imshow(uint8(img2));
title('Image #2','FontSize',16);hold on;
for n=1:size(xCorner2,1)
    rectangle('Position',[xCorner2(n)-4,yCorner2(n)-4,9,9],'LineWidth',1.1,'EdgeColor','y');
end
hold off;

%% (b) Feature matching %%
winSize_corr=9;p=floor(winSize_corr/2);
Thres_corr=0.84;
rmBorder1=find((round(xCorner1)-p>1)&(round(yCorner1)-p>1)&(round(xCorner1)+p<size(img1,2))&(round(yCorner1)+p<size(img1,1)));
rmBorder2=find((round(xCorner2)-p>1)&(round(yCorner2)-p>1)&(round(xCorner2)+p<size(img2,2))&(round(yCorner2)+p<size(img2,1)));
xCorner1_rb=xCorner1(rmBorder1);yCorner1_rb=yCorner1(rmBorder1);
xCorner2_rb=xCorner2(rmBorder2);yCorner2_rb=yCorner2(rmBorder2);

%[xMatch1,yMatch1,xMatch12,yMatch12]=nCorrCoef(img1,img2,xCorner1_rb,yCorner1_rb,xCorner2_rb,yCorner2_rb,Thres_corr,winSize_corr);
%save('nCorrCoef_new.mat','xMatch1','yMatch1','xMatch12','yMatch12');
%load('nCorrCoef_new');

%[xMatch1,yMatch1,xMatch12,yMatch12]=nCorrCoef_interp(img1,img2,xCorner1_rb,yCorner1_rb,xCorner2_rb,yCorner2_rb,Thres_corr,winSize_corr);
%save('nCorrCoef_interp.mat','xMatch1','yMatch1','xMatch12','yMatch12');
load('nCorrCoef_interp');

figure(3);
subplot(121);imshow(uint8(img1));
title('Image #1','FontSize',16);hold on;
for n1=1:numel(xMatch1)
    rectangle('Position',[xMatch1(n1)-4,yMatch1(n1)-4,9,9],'LineWidth',1.1,'EdgeColor','y');
    plot([xMatch1(n1);xMatch12(n1)],[yMatch1(n1);yMatch12(n1)],'y','LineWidth',1.1);
end
hold off;
subplot(122);imshow(uint8(img2));
title('Image #2','FontSize',16);hold on;
for n2=1:numel(xMatch12)
    rectangle('Position',[xMatch12(n2)-4,yMatch12(n2)-4,9,9],'LineWidth',1.1,'EdgeColor','y');
    plot([xMatch12(n2);xMatch1(n2)],[yMatch12(n2);yMatch1(n2)],'y','LineWidth',1.1);
end
hold off;

%% (c) Outlier rejection (MSAC: 7-point algorithm) %%
%[xIn1,yIn1,xIn12,yIn12,numInliers,numTrials,dltIn]=MSAC(xMatch1,yMatch1,xMatch12,yMatch12);   
load('MSACtest_interp2');

figure(4);
subplot(121);imshow(uint8(img1));
title('Image #1','FontSize',16);hold on;
for n1=1:numel(xIn1)
    rectangle('Position',[xIn1(n1)-4,yIn1(n1)-4,9,9],'LineWidth',1.1,'EdgeColor','y');
    plot([xIn1(n1);xIn12(n1)],[yIn1(n1);yIn12(n1)],'y','LineWidth',1.1);
end
hold off;
subplot(122);imshow(uint8(img2));
title('Image #2','FontSize',16);hold on;
for n2=1:numel(xIn12)
    rectangle('Position',[xIn12(n2)-4,yIn12(n2)-4,9,9],'LineWidth',1.1,'EdgeColor','y');
    plot([xIn12(n2);xIn1(n2)],[yIn12(n2);yIn1(n2)],'y','LineWidth',1.1);
end
hold off;

toc;

%% (d) Linear estimation (DLT) %%
x1=[xIn1,yIn1];x12=[xIn12,yIn12];
[F_DLT]=DLT(x1,x12,numInliers);

%% (e) Nonlinear estimation (L-M) %%
% F_DLT=F_DLT';F_DLT=F_DLT(:);F_DLT=Paramtrz(F_DLT);F_DLT=deParamtrz(F_DLT);F_DLT=reshape(F_DLT,3,3)';
% F_DLT=-F_DLT;
%% STEP #1. Initialization
P=[eye(3),zeros(3,1)];
Pp=getPp(F_DLT);
Pp=Pp/norm(Pp,'fro');

X_pi=Triangulation( xIn1,yIn1,xIn12,yIn12,F_DLT,Pp ); % initialize 3D points

% Parameter Vector:
Pp_T=Pp';php=Paramtrz(Pp_T(:));
Xh=zeros(3,numInliers);
for pt=1:numInliers
    Xh(:,pt)=Paramtrz(X_pi(:,pt));
end
Xh=Xh(:);
PV=cat(2,php',Xh')';

% Measurement Vector:
x=[xIn1';yIn1'];
xp=[xIn12';yIn12'];
MV=cat(2,x(:)',xp(:)')';

CST=[];
iterationTimes=0;
n=numInliers;
lbd=0.001;

xh=homo2inhomo(P*X_pi);xhp=homo2inhomo(Pp*X_pi);

EP=x-xh;EPP=xp-xhp;
%ep = vet(EP) (2n x 1):
ep=zeros(2*n,1);epp=zeros(2*n,1);
ep(1:2:2*n-1)=EP(1,:);epp(1:2:2*n-1)=EPP(1,:);
ep(2:2:2*n)=EP(2,:);epp(2:2:2*n)=EPP(2,:);

SIGMAx=eye(2*n);
cst=ep'*inv(SIGMAx)*ep+epp'*inv(SIGMAx)*epp;
CST(1)=cst;
%% STEP #2. Jacobian Matrix
%while 1
for it=1:5    
    [~,Pp,Xh,~]=getFrPV(PV);

    Aip=cell(1,numInliers);Bi=cell(1,numInliers);Bip=cell(1,numInliers);
    for i=1:numInliers
        %Aip{i}=jacobianA(xhp_hm(:,i),php);
        Aip{i}=jacobianMatrixA(Xh(:,i),Pp);
        Bi{i}=jacobianB(Xh(:,i),Pp);
        Bip{i}=jacobianB(Xh(:,i),P);
    end
    %% STEP #3. Normal Equation
    Up=zeros(11);SIGMAxip=eye(2);
    V=cell(1,n);
    Wp=cell(1,n);
    for i=1:n
        Up=Up+Aip{i}'*inv(SIGMAxip)*Aip{i}; % 11x11
        V{i}=Bi{i}'*inv(SIGMAxip)*Bi{i}+Bip{i}'*inv(SIGMAxip)*Bip{i}; % 3x3
        Wp{i}=Aip{i}'*inv(SIGMAxip)*Bip{i}; % 11x3
    end

        %Normal Equation Vector
        ep_ap=zeros(11,1);ep_b=zeros(3,n);
        [~,Pp,~,X_pi]=getFrPV(PV);
        xh=homo2inhomo(P*X_pi);xhp=homo2inhomo(Pp*X_pi);
        EP=x-xh;EPP=xp-xhp;
        for i=1:n
            ep_ap=ep_ap+Aip{i}'*inv(SIGMAxip)*EPP(:,i); % 11x1
            ep_b(:,i)=Bi{i}'*inv(SIGMAxip)*EP(:,i) + Bip{i}'*inv(SIGMAxip)*EPP(:,i); % 3x1
        end
    %% STEP #4 update delta
    while 1
       Up_star=Up+lbd*eye(size(Up)); % augmentated Upp
       V_star=cell(size(V));
       Sp=Up_star;ep_Big=ep_ap;
       for i=1:n
           V_star{i}=V{i}+lbd*eye(size(V{i}));
           Sp=Sp-Wp{i}*inv(V_star{i})*Wp{i}';
           ep_Big=ep_Big-Wp{i}*inv(V_star{i})*ep_b(:,i);
       end
       delta_ap=Sp\ep_Big;
       delta_b=zeros(3,n);
       for i=1:n
           delta_b(:,i)=V_star{i}\(ep_b(:,i)-Wp{i}'*delta_ap);
       end
       delta=[delta_ap;delta_b(:)];
       %% STEP #5. Candidate Parameter Vector
       PV0=PV+delta;
       %% STEP #6. update errors
       [~,Pp0,~,X_pi0]=getFrPV(PV0);
       xh0=homo2inhomo(P*X_pi0);xh0p=homo2inhomo(Pp0*X_pi0);
       EP0=x-xh0;EP0P=xp-xh0p;
       e0p=zeros(2*n,1);e0pp=zeros(2*n,1);
       e0p(1:2:2*n-1)=EP0(1,:);e0pp(1:2:2*n-1)=EP0P(1,:);
       e0p(2:2:2*n)=EP0(2,:);e0pp(2:2:2*n)=EP0P(2,:);
       %% STEP #7. update cost and determine where to go next 
       cst0=e0p'*inv(SIGMAx)*e0p+e0pp'*inv(SIGMAx)*e0pp;
       if cst0 < cst
           PV=PV0;
           cst=cst0;
           ep=e0p;epp=e0pp;
           lbd=0.1*lbd;

           CST=[CST,cst0];
           iterationTimes=iterationTimes+1;
           break;
       else
           lbd=10.0*lbd;
       end
    end
end
%% get F_LM
[~,Pp_LM,~,~]=getFrPV(PV);
M=Pp_LM(1:3,1:3);Ep=Pp_LM(:,end);
F_LM=skewSymMatrix(Ep)*M;
F_LM=F_LM/norm(F_LM,'fro')

CST'
figure(5); plot(0:iterationTimes,CST,'b-x','LineWidth',0.7);grid on;
xlabel('iterations','FontSize',15);ylabel('costs','FontSize',15);
title('The Trend of Cost vs. Iteration for the Levenberg-Marquardt Algorithm','FontSize',15);

%% (f) Point to line mapping %%
pk3Pts=pickOutliers( xMatch1,yMatch1,xIn1,yIn1 );
%pk3Pts=[66,67,158];
x_3=[xMatch1(pk3Pts),yMatch1(pk3Pts)]';

figure(6);
subplot(121);imshow(uint8(img1));
title('Image #1','FontSize',16);hold on;

for pt=1:3
    rectangle('Position',[x_3(1,pt)-4,x_3(2,pt)-4,9,9],'LineWidth',1.1,'EdgeColor','r');
end
hold off;

lp_3=F_LM*inhomo2homo(x_3);

subplot(122);imshow(uint8(img2));
title('Image #2','FontSize',16);hold on;

for li=1:3
    cefs=[-lp_3(1,li)/lp_3(2,li),-lp_3(3,li)/lp_3(2,li)];
    h=refline(cefs);
    set(h,'Color','r');
end

hold off;
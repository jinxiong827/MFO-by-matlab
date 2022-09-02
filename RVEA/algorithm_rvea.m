clear;clc;close all
alpha=2;
mu=5;
wmax=100;
Problem.M=2;
Problem.D=20;
Problem.N=100;
maxFE=300;
FE=0;
Problem.name='DTLZ2';
Problem.upper=ones(1,Problem.D);
Problem.lower=zeros(1,Problem.D);
[V0,Problem.N] = UniformPoint2(Problem.N,Problem.M);
V     = V0;
NI    = 11*Problem.D-1;
P     = UniformPoint2(NI,Problem.D,'Latin');
[A2,FE]    = SOLUTION2(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1),FE,Problem.M,Problem.name);
A1    = A2;
THETA = 5.*ones(Problem.M,Problem.D);
Model = cell(1,Problem.M);

while FE<maxFE
    A1Dec = cat(1,A1.dec);
    A1Obj = cat(1,A1.obj);

%     A1Dec = unique(A1Dec,'rows');
%     A1Obj = unique(A1Obj,'rows');
    for i = 1 : Problem.M
        dmodel     = dacefit(A1Dec,A1Obj(:,i),'regpoly1','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
        Model{i}   = dmodel;
        THETA(i,:) = dmodel.theta;
    end
    PopDec = A1Dec;
    w      = 1;

    while w <= wmax
        drawnow();
        OffDec = OperatorGA2(PopDec);
        PopDec = [PopDec;OffDec];
        [N,~]  = size(PopDec);
        PopObj = zeros(N,Problem.M);
        MSE    = zeros(N,Problem.M);
        for i = 1: N
            for j = 1 : Problem.M
                [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),Model{j});
            end
        end
        index  = KEnvironmentalSelection2(PopObj,V,(w/wmax)^alpha);
        PopDec = PopDec(index,:);
        PopObj = PopObj(index,:);
        % Adapt referece vectors
        if ~mod(w,ceil(wmax*0.1))
            V(1:Problem.N,:) = V0.*repmat(max(PopObj,[],1)-min(PopObj,[],1),size(V0,1),1);
        end
        w = w + 1;
        %                 plot(PopObj(:,1),PopObj(:,2),'*')
    end

    [NumVf,~] = NoActive2(A1Obj,V0);
    PopNew    = KrigingSelect2(PopDec,PopObj,MSE(index,:),V,V0,NumVf,0.05*Problem.N,mu,(w/wmax)^alpha);
    [New,FE]  = SOLUTION2(PopNew,FE,Problem.M,Problem.name);
    A2        = [A2,New];
    A1        = UpdataArchive2(A1,New,V,mu,NI);

    A2Obj=cat(1,A2.obj);
    if Problem.M==2
    plot(A2Obj(:,1),A2Obj(:,2),'*')
    elseif Problem.M==3
    plot3(A2Obj(:,1),A2Obj(:,2),A2Obj(:,3),'*')    
    end
end
A2Best=A2Obj(NDSort2(cat(1,A2.obj),1)==1,:);
pf=PF(Problem,10000);
pf_2=PF(Problem,1000);
hold on
if Problem.M==2
plot(pf(:,1),pf(:,2),'-')
elseif Problem.M==3
plot3(pf_2(:,1),pf_2(:,2),pf_2(:,3),'*')
end
Final_IGD=IGD2(A2Best,pf);
Final_IGD
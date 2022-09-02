clear;clc;close all
Ke=5;
delta=0.9;
nr=2;
L1=80;
L2=20;
maxFE=110;
FE=0;
Problem.M=3;
Problem.D=10;
Problem.N=100;
Problem.name='DTLZ7';
Problem.upper=ones(1,Problem.D);
Problem.lower=zeros(1,Problem.D);

NI = 11*Problem.D-1;
P  = UniformPoint2(NI,Problem.D,'Latin');
[Population,FE]    = SOLUTION2(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1),FE,Problem.M,Problem.name);
L1 = min(L1,length(Population));

while FE<maxFE
    drawnow();
    PopDec    = EGOSelect2(Problem,Population,L1,L2,Ke,delta,nr);
    [Offspring,FE]  = SOLUTION2(PopDec,FE,Problem.M,Problem.name);
    Population = [Population,Offspring];
    PopulationObj=cat(1,Population.obj);

    if Problem.M==2
    plot(PopulationObj(:,1),PopulationObj(:,2),'*')
    elseif Problem.M==3
    plot3(PopulationObj(:,1),PopulationObj(:,2),PopulationObj(:,3),'*')   
    end
end

PopulationObjBest=PopulationObj(NDSort2(cat(1,Population.obj),1)==1,:);
pf=PF(Problem,10000);
pf_2=PF(Problem,1000);
hold on
if Problem.M==2
plot(pf(:,1),pf(:,2),'-')
elseif Problem.M==3
plot3(pf_2(:,1),pf_2(:,2),pf_2(:,3),'*')
end
Final_IGD=IGD2(PopulationObjBest,pf);
Final_IGD

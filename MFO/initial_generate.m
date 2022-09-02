clc;clear

%%
No_initial=150;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
No_blade=21;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
No_blade_mid=11;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
No_blade_point=100;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
profile='NACA';%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ro=1.2;
No_fan=10;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

Pt_design=[25	100	25	100	25	100	25 50 150	150]; %static pressure (pa)%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Q_design=[40000	60000	15000	70000	20000	80000	10000	40000	90000	150000]; %volumn flow (m3/h)%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
D2_design=[1	1	0.8	1.1	0.9	1.2	0.7	1.29	1.3	1.6]; %diameter (m)%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
N_design=[900	1400	750	1400	750	1400	900	600	1400	1400];%rotation speed%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

blade_filepath='/test4T/xiongjin/57fanopt/56';% running path
pre_filepath='\ansys_inc\v190\CFX\bin'; %ansys path
pre_filepathl=replace(pre_filepath,'\','/');
Zin=200;
Zout=-200;
tolerance=[2	2	2	2.5	2.5	2.5	3	3	3	4];

%% scale
if size(Q_design,2)~=No_fan || size(D2_design,2)~=No_fan || size(N_design,2)~=No_fan
    error('wrong input')
end
u2=pi.*D2_design.*N_design/60;
r2=D2_design/2;
% Pt_design=P_design+0.5*(Q_design*4/3600./(pi.*D2_design.^2)).^2*ro;
Qbar=Q_design./3600./((pi/4).*D2_design.^2.*u2);%flow coefficient
Pbar=Pt_design/ro./u2.^2;%pressure coefficient
Nbar=1000*N_design./((pi/4).*D2_design.^2*ro.*u2.^5);
ns=14.8*Qbar.^0.5./Pbar.^0.75;
Ns=N_design.*(Q_design./3600).^0.5./Pt_design.^0.75;
%  7.8137
N_median=median_value(N_design);
D2_median=median_value(D2_design);
tolerance_median=median_value(tolerance);
u2_median=pi*D2_median*N_median/60;
Q_scale=Qbar*(pi/4)*u2_median*D2_median^2;
Q_scale=Q_scale*1.03
% Qm_scale=Q_scale*ro;
P_scale=Pbar*ro*u2_median^2
criterion=sortrows([P_scale;Q_scale]',-2);
%  criterion=[P_scale;Q_scale]';
criterion=[1:No_fan;criterion']'
%search zone
Qmax=criterion(1,3)*1.02;
Qmin=criterion(end,3)*0.98;
Pmax=max(criterion(:,2))*1.02;
Pmin=min(criterion(:,2))*0.98;
%
plot(criterion(:,2),criterion(:,3),'*')
A_text={'1','2','3','4','5','6','7','8','9','10'};
text(criterion(:,2)+0.001,criterion(:,3)-0.001,A_text);
hold on
%%
%initial_design & parametrization
p=11; %
ABC=lhsdesign(No_initial,p,'criterion','correlation')';
% ABC=lhsdesign(No_initial,p)';
%C1
c1min=Qmin;c1max=Qmax;
c1=ABC(1,:)*(c1max-c1min)+c1min;%C1:Qin=[Qmin,Qmax]
%C2
c2min=0.1;c2max=0.3;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c2=ABC(2,:)*(c2max-c2min)+c2min;%C2: D1/D2=[0.2~0.4]
D1=D2_median*c2;
for i=1:length(D1)
    D_No_blade(i,:)=linspace(D1(i),D2_median,No_blade);%diameter for each blade section
    r_No_blade=D_No_blade/2;
end
D1_blade=D1*0.9;
%C3~C5
%Angel_incidence
cz=c1./(pi/4*(D2_median^2-D1.^2));
u2_shroud=u2_median;
u2_mid=u2_median/D2_median*(D1+D2_median)/2;
u2_hub=u2_median/D2_median*D1;
Angel_incidence=4;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c3min=0;c3max=2*Angel_incidence;
c3=ABC(3,:)*(c3max-c3min)+c3min;%C3: i_shroud
c4min=-Angel_incidence;c4max=Angel_incidence;
c4=ABC(4,:)*(c4max-c4min)+c4min;%C4: i_mid
c5min=-3*Angel_incidence;c5max=-Angel_incidence;
c5=ABC(5,:)*(c5max-c5min)+c5min;%C5: i_hub

beta_1_shroud=atand(cz./u2_shroud)+c3;
beta_1_mid=atand(cz./u2_mid)+c4;
beta_1_hub=atand(cz./u2_hub)+c5;
for i=1:size(ABC,2)
    beta_1A_c(i,:)=polyfit([1 No_blade_mid No_blade],[beta_1_hub(i) beta_1_mid(i) beta_1_shroud(i)],2);
    beta_1A(i,1:No_blade)=polyval(beta_1A_c(i,:),1:No_blade);% beta1A hub->shroud
end
% beta_1A
%C6~C7
%C6 total pressure rise
%C7 coefficient of free-vortex fan
% beta_2profile

c6min=Pmin;c6max=Pmax;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c6=ABC(6,:)*(c6max-c6min)+c6min;%C6=[Pmin Pmax]: total pressure rise
P_theoretical=c6/0.8;
c7min=0;c7max=1;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c7=ABC(7,:)*(c7max-c7min)+c7min;%C7=[-1 1]: coefficient of free-vortex fan

D_geom_mean=((D1.^2+D2_median^2)/2).^0.5;
% D_geom_mean=(D1+D2_median)/2;
u2_geom_mean=u2_median/D2_median*D_geom_mean;
dCu_m=P_theoretical./u2_geom_mean/ro;
rdCu_m=dCu_m.*D_geom_mean.^c7;
for i=1:length(rdCu_m)
    dCu(i,:)=rdCu_m(i)./D_No_blade(i,:).^c7(i); %dCu hub->shroud
end

u2_No_blade=D_No_blade*pi*N_median/60;
P_theoretical_No_blade=ro*u2_No_blade.*dCu;
cz_No_blade=u2_No_blade.*tand(beta_1A);
cz_ave=sum(cz_No_blade.*D_No_blade.^2,2)./sum(D_No_blade.^2,2);
cz_ave=repmat(cz_ave,1,No_blade);
beta_2=atand(cz_ave./(u2_No_blade-dCu));
beta_2(beta_2<0)=beta_2(beta_2<0)+180;
behind_angle=3;%<<<<<<<<<<<<<<<<<<behind_angle=3
beta_2A=beta_2+behind_angle;

%C8
c8min=0.8;c8max=1.3;%<<<<<<<<<<<<<<<<<<c8:cy_est=[0.8 1.8]
c8=ABC(8,:)*(c8max-c8min)+c8min;%
cy=repmat(c8',1,No_blade);
% t_No_blade=pi*D_No_blade/Z;
%calculate b1
wm_No_blade=(cz_No_blade.^2+(u2_No_blade-dCu).^2).^0.5;
% b_1_No_blade=2*t_No_blade.*P_theoretical_No_blade./(ro.*u2_No_blade.*wm_No_blade.*cy)
c9min=3;c9max=6.999;%<<<<<<<<<<<<<<<<<<
c9=ABC(9,:)*(c9max-c9min)+c9min;%
Z=fix(repmat(c9',1,No_blade));%blade number<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
b_No_blade=4*pi.*P_theoretical_No_blade./(cy.*Z.*ro.*(2*pi*N_median/60).*wm_No_blade);
b_No_blade=b_No_blade*1.1;
for i=1:fix(No_blade/4)
    b_No_blade(:,i)=b_No_blade(:,No_blade-(i-1)*3);
end
 b_No_blade(:,6)=1.5*b_No_blade(:,5)-0.5*b_No_blade(:,4);
 b_No_blade(:,8)=1.5*b_No_blade(:,9)-0.5*b_No_blade(:,10);
 b_No_blade(:,7)=(1.5*b_No_blade(:,6)-0.5*b_No_blade(:,5)+1.5*b_No_blade(:,9)-0.5*b_No_blade(:,10))/2;
%C9~C11
c10min=0.3;c10max=0.8;%<<<<<<<<<<<<<<<<<<c9:a/b=[0.3 0.6]
c10=ABC(9,:)*(c10max-c10min)+c10min;
c11min=0.06;c11max=0.24;%<<<<<<<<<<<<<<<<<<max blade thickness:[6% 24%]
c11=ABC(10,:)*(c11max-c11min)+c11min;
max_t=c11; % max blade thickness

% [~,maxb_index]=max(b_No_blade')
% for ii=1:50
%     plot(1:No_blade,b_No_blade(ii,:))
%     hold on
% end

Ge=1;
diary on
flag0=zeros(1,No_initial);
for i=1:size(ABC,2)
    file1=sprintf('%s%s%d%s',blade_filepath,'/blade',Ge*1000+i,'.crv');
    fid1=fopen(file1,'w');
    fprintf(fid1,'%s\n','## ExportPoints1');
    for j=1:No_blade
        meanline_xy=meanline_parabola(b_No_blade(i,j),beta_1A(i,j),beta_2(i,j),c10(i),No_blade_point);
        blade_profile_xy=blade_profile(meanline_xy,b_No_blade(i,j),max_t(i),profile)*1000;
        fprintf(fid1,'%s%d%s\n','# Profile 1 at    ',(j-1)*100/(No_blade-1),'.0000%');
        fprintf(fid1,'%.8f\t%.8f\t%.8f   le1\n',blade_profile_xy(1,1),((D2_median/2-D1_blade(i)/2)/(No_blade-1)*(j-1)+D1_blade(i)/2)*1000,blade_profile_xy(2,1));
        for k=2:size(blade_profile_xy,2)
            fprintf(fid1,'%.8f\t%.8f\t%.8f \n',blade_profile_xy(1,k),((D2_median/2-D1_blade(i)/2)/(No_blade-1)*(j-1)+D1_blade(i)/2)*1000,blade_profile_xy(2,k));
        end
        fprintf(fid1,'%.8f\t%.8f\t%.8f   te1\n',blade_profile_xy(1,end)/2+blade_profile_xy(3,end)/2,((D2_median/2-D1_blade(i)/2)/(No_blade-1)*(j-1)+D1_blade(i)/2)*1000,blade_profile_xy(2,end)/2+blade_profile_xy(4,end)/2);
        for k=1:size(blade_profile_xy,2)-1
            fprintf(fid1,'%.8f\t%.8f\t%.8f \n',blade_profile_xy(3,end-k+1),((D2_median/2-D1_blade(i)/2)/(No_blade-1)*(j-1)+D1_blade(i)/2)*1000,blade_profile_xy(4,end-k+1));
        end
        
        if flag0(i)<max(abs([blade_profile_xy(2,:) blade_profile_xy(4,:)]))
            flag0(i)=max(abs([blade_profile_xy(2,:) blade_profile_xy(4,:)]));
        end
    end
    fclose(fid1);
    if flag0(i)<=180
        
        file1=sprintf('%s%s%d%s',blade_filepath,'/shroud',Ge*1000+i,'.crv');
        fid1=fopen(file1,'w');
        for j=1:500
            fprintf(fid1,'%.8f\t%.8f\t%.8f\r\n',0,D2_median/2*1000-0.1,(Zout-Zin)/499*(j-1)+Zin);
        end
        fclose(fid1);
        
        file1=sprintf('%s%s%d%s',blade_filepath,'/hub',Ge*1000+i,'.crv');
        fid1=fopen(file1,'w');
        for j=1:500
            fprintf(fid1,'%.8f\t%.8f\t%.8f\r\n',0,D1(i)/2*1000,(Zout-Zin)/499*(j-1)+Zin);
        end
        fclose(fid1);
        
        turbogrid(Ge,i,tolerance_median,600000,blade_filepath,blade_filepath);
        fid0=fopen('cfx_turbogrid.sh','w');
        fprintf(fid0,'%s%d%s%d%s\n','/ansys_inc/v190/TurboGrid/bin/cfxtg -batch /test4T/xiongjin/57fanopt/',56,'/blade',1000*Ge+i,'.tse');
        fclose(fid0);
        fprintf('%s%d\n','now running: ',Ge*1000+i);
                !sh cfx_turbogrid.sh
        fid4=fopen('cfxpre.sh','w');
        cfxpre(Ge,i,N_median,c1(i)/Z(i,1),blade_filepath,blade_filepath)
        fprintf(fid4,'%s%d%s%d%s\n','/ansys_inc/v190/CFX/bin/cfx5pre -batch /test4T/xiongjin/57fanopt/',56,'/blade',1000*Ge+i,'.pre');
        fclose(fid4);
        pause(5);
                !sh cfxpre.sh
        pause(5);
                !pkill SolverManager
    end
end
diary off
%==========================
fid1=fopen('cfx_solve.sh','w');
for i=1:size(ABC,2)
    fprintf(fid1,'%s%d%s\n','echo "',1000*Ge+i,' start"');
    fprintf(fid1,'%s%d%s%d%s\n','/ansys_inc/v190/CFX/bin/cfx5solve -def /test4T/xiongjin/57fanopt/',56,'/blade',1000*Ge+i,'.def -part 115 -start-method "Platform MPI Local Parallel" ;');
end
fclose(fid1);
save data01
!sh cfx_solve.sh
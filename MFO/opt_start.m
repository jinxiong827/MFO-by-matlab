clc;clear;close all
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
tolerance=[2	2	2	2.5	2.5	2.5	3	3	3	4];% tip tolerance distance%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
%
N_median=median_value(N_design);
D2_median=median_value(D2_design);
tolerance_median=median_value(tolerance);
u2_median=pi*D2_median*N_median/60;
Q_scale=Qbar*(pi/4)*u2_median*D2_median^2;
Q_scale=Q_scale*1.03;
% Qm_scale=Q_scale*ro;
P_scale=Pbar*ro*u2_median^2;
criterion=sortrows([P_scale;Q_scale]',-2);
criterion=[1:No_fan;criterion']';
%search zone
Qmax=criterion(1,3)*1.02;
Qmin=criterion(end,3)*0.98;
Pmax=max(criterion(:,2))*1.02;
Pmin=min(criterion(:,2))*0.98;
%==================================foregoing code is same as initial_generate.m
initial_design=[]; %data from initial_generate.m 

criterion=[1   25.0000   11.4444
    2   41.3265   11.0357
    3   56.2500   10.0586
    4   34.1542    9.6732
    5   44.4444    9.4193
    6   28.6990    8.5152
    7   51.0204    8.3414
    8   36.6804    7.5346
    9   24.2148    6.7357
    10   69.4444    6.6229];


p=11;

new_load=load('new_result.txt');
new_load0(:,1:p)=new_load(:,2:p+1);
new_load0(:,p+1:p+2)=new_load(:,17:18);
new_load0(:,p+3)=new_load(:,13);
new_load0(:,p+4)=new_load(:,25);
initial_design=[initial_design;new_load0];

Q_c1_design=(criterion(:,3)-Qmin)/(Qmax-Qmin);

father=initial_design;
father_design=initial_design(:,p);
father_pressure=initial_design(:,p+1);
father_efficiency=initial_design(:,p+2);
father_static_pressure=initial_design(:,p+3);
father_static_efficiency=initial_design(:,p+4);

Afmax=dir(strcat(blade_filepath,'/blade*.res'));
for i=1:size(Afmax,1)
    res_file_max(i)=str2num(Afmax(i).name(6:9));
end
No=size(criterion,1);
% Ge=2;
% kn=1;
Ge=fix(max(res_file_max)/1000);
kn=(fix((max(res_file_max)-Ge*1000)/No)+1)*No;
flag1=0;
theta = 0.5*ones(1,p); lob = 0.00001*ones(1,p); upb = ones(1,p);
min_range=0*ones(1,p-1);
max_range=ones(1,p-1);
while(1)
    
    if kn<999-No
        kn=kn+No;
    else
        Ge=Ge+1;
        kn=1;
    end
    
    
    [mS,mY]=dsmerge(father(:,1:p),father(:,p+3:p+4));
    if flag1==0
        [dmodel01, ~] =dacefit(mS, mY(:,1), @regpoly0, @corrgauss, theta, lob, upb);
        [dmodel02, ~] =dacefit(mS, mY(:,2), @regpoly0, @corrgauss, theta, lob, upb);
    elseif flag1==1
        [dmodel01, ~] =dacefit(mS, mY(:,1), @regpoly1, @corrgauss, theta, lob, upb);
        [dmodel02, ~] =dacefit(mS, mY(:,2), @regpoly1, @corrgauss, theta, lob, upb);
    elseif  flag1==2
        [dmodel01, ~] =dacefit(mS, mY(:,1), @regpoly2, @corrgauss, theta, lob, upb);
        [dmodel02, ~] =dacefit(mS, mY(:,2), @regpoly2, @corrgauss, theta, lob, upb);
    end
    save modeldata01 dmodel01
    save modeldata02 dmodel02
    x00=[];
    fval=[];
    for i=1:size(criterion,1)
        
        model_Q=Q_c1_design(i);
        model_P=criterion(i,2);
        save modeldata_c  model_Q model_P
        fitnessfcn=@object;
        options=gaoptimset('paretoFraction',30/100,'populationsize',100,'generations',300,'stallGenLimit',200,'TolFun',1e-100);
        [x00,fval]=gamultiobj(fitnessfcn,p-1,[],[],[],[],min_range,max_range,options);
        y00=sortrows([x00 fval],p);
        y00=y00(1:sum(y00(:,p)<model_P*0.05),:);
        new_offspring(i,1)=model_Q;new_offspring_b(i,1)=model_Q;new_offspring_c(i,1)=model_Q;
        y00_c=y00(find(y00(:,p+1)==min(y00(:,p+1))),1:p-1);
        y00_b=y00(1,1:p-1);
        
        if max(y00(2,1:p-1)-y00_b)<0.01 || max(y00(2,1:p-1)-y00_c)<0.01
            y00_d=rand(1,p-1);
        else
            y00_d=y00(2,1:p-1);
        end
        
        dis_new_father=10000;
        for j=1:size(father,1)
            if max(father(j,2:11)-y00_b)<dis_new_father
                dis_new_father=max(father(j,2:11)-y00_b);
            end
        end
        
        if dis_new_father<0.05
            new_offspring(i,2:p)=y00_c(1,:);
            new_offspring_b(i,2:p)=y00_d(1,:);
            new_offspring_c(i,2:p)=rand(1,p-1);
        else
            new_offspring(i,2:p)=y00_b(1,:);
            new_offspring_b(i,2:p)=y00_c(1,:);
            new_offspring_c(i,2:p)=y00_d(1,:);
        end
        
       
        for i2=1:4
            
            if i2==1
                ABC(:,1)=new_offspring(i,:)';
            elseif i2==2
                ABC(:,1)=new_offspring_b(i,:)';
            elseif i2==3
                ABC(:,1)=new_offspring_c(i,:)';
            else
                ABC(1,1)=new_offspring(i,1);
                ABC(2:p,1)=rand(p-1,1);
            end
            
            % ABC=lhsdesign(No_initial,p)';
            %C1
            c1min=Qmin;c1max=Qmax;
            c1=ABC(1,:)*(c1max-c1min)+c1min;%C1:Qin=[Qmin,Qmax]
            %C2
            c2min=0.1;c2max=0.3;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            c2=ABC(2,:)*(c2max-c2min)+c2min;%C2: D1/D2=[0.2~0.4]
            D1=D2_median*c2;
            for i1=1:length(D1)
                D_No_blade(i1,:)=linspace(D1(i1),D2_median,No_blade);%diameter for each blade section
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
            for i1=1:size(ABC,2)
                beta_1A_c(i1,:)=polyfit([1 No_blade_mid No_blade],[beta_1_hub(i1) beta_1_mid(i1) beta_1_shroud(i1)],2);
                beta_1A(i1,1:No_blade)=polyval(beta_1A_c(i1,:),1:No_blade);% beta1A hub->shroud
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
            for i1=1:length(rdCu_m)
                dCu(i1,:)=rdCu_m(i1)./D_No_blade(i1,:).^c7(i1); %dCu hub->shroud
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
            c9min=3;c9max=6.999;%<<<<<<<<<<<<<<<<<<c8:cy_est=[0.8 1.8]
            c9=ABC(9,:)*(c9max-c9min)+c9min;%
            Z=fix(repmat(c9',1,No_blade));%blade number<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            b_No_blade=4*pi.*P_theoretical_No_blade./(cy.*Z.*ro.*(2*pi*N_median/60).*wm_No_blade);
            b_No_blade=b_No_blade*1.1;
            for i1=1:fix(No_blade/4)
                b_No_blade(:,i1)=b_No_blade(:,No_blade-(i1-1)*3);
            end
            b_No_blade(:,6)=1.5*b_No_blade(:,5)-0.5*b_No_blade(:,4);
            b_No_blade(:,8)=1.5*b_No_blade(:,9)-0.5*b_No_blade(:,10);
            b_No_blade(:,7)=(1.5*b_No_blade(:,6)-0.5*b_No_blade(:,5)+1.5*b_No_blade(:,9)-0.5*b_No_blade(:,10))/2;
            
            %C9~C11
            c10min=0.3;c10max=0.8;%<<<<<<<<<<<<<<<<<<c9:a/b=[0.3 0.6]
            c10=ABC(9,:)*(c10max-c10min)+c10min;
            c11min=0.06;c11max=0.18;%<<<<<<<<<<<<<<<<<<max blade thickness:[6% 24%]
            c11=ABC(10,:)*(c11max-c11min)+c11min;
            max_t=c11; % max blade thickness
            
            
            flag0=zeros(1,No_blade);
            file1=sprintf('%s%s%d%s',blade_filepath,'/blade',Ge*1000+i+kn-1,'.crv');
            fid1=fopen(file1,'w');
            fprintf(fid1,'%s\n','## ExportPoints1');
            for j=1:No_blade
                meanline_xy=meanline_parabola(b_No_blade(1,j),beta_1A(1,j),beta_2(1,j),c10(1),No_blade_point);
                blade_profile_xy=blade_profile(meanline_xy,b_No_blade(1,j),max_t(1),profile)*1000;
                fprintf(fid1,'%s%d%s\n','# Profile 1 at    ',(j-1)*100/(No_blade-1),'.0000%');
                fprintf(fid1,'%.8f\t%.8f\t%.8f   le1\n',blade_profile_xy(1,1),((D2_median/2-D1_blade(1)/2)/(No_blade-1)*(j-1)+D1_blade(1)/2)*1000,blade_profile_xy(2,1));
                for k=2:size(blade_profile_xy,2)
                    fprintf(fid1,'%.8f\t%.8f\t%.8f \n',blade_profile_xy(1,k),((D2_median/2-D1_blade(1)/2)/(No_blade-1)*(j-1)+D1_blade(1)/2)*1000,blade_profile_xy(2,k));
                end
                fprintf(fid1,'%.8f\t%.8f\t%.8f   te1\n',blade_profile_xy(1,end)/2+blade_profile_xy(3,end)/2,((D2_median/2-D1_blade(1)/2)/(No_blade-1)*(j-1)+D1_blade(1)/2)*1000,blade_profile_xy(2,end)/2+blade_profile_xy(4,end)/2);
                for k=1:size(blade_profile_xy,2)-1
                    fprintf(fid1,'%.8f\t%.8f\t%.8f \n',blade_profile_xy(3,end-k+1),((D2_median/2-D1_blade(1)/2)/(No_blade-1)*(j-1)+D1_blade(1)/2)*1000,blade_profile_xy(4,end-k+1));
                end
                
                if flag0(j)<max(abs([blade_profile_xy(2,:) blade_profile_xy(4,:)]))
                    flag0(j)=max(abs([blade_profile_xy(2,:) blade_profile_xy(4,:)]));
                end
            end
            fclose(fid1);
            if max(flag0)<=180
                file1=sprintf('%s%s%d%s',blade_filepath,'/shroud',Ge*1000+i+kn-1,'.crv');
                fid1=fopen(file1,'w');
                for j=1:500
                    fprintf(fid1,'%.8f\t%.8f\t%.8f\r\n',0,D2_median/2*1000-0.1,(Zout-Zin)/499*(j-1)+Zin);
                end
                fclose(fid1);
                file1=sprintf('%s%s%d%s',blade_filepath,'/hub',Ge*1000+i+kn-1,'.crv');
                fid1=fopen(file1,'w');
                for j=1:500
                    fprintf(fid1,'%.8f\t%.8f\t%.8f\r\n',0,D1(i)/2*1000,(Zout-Zin)/499*(j-1)+Zin);
                end
                fclose(fid1);
                turbogrid(Ge,i+kn-1,tolerance_median,600000,blade_filepath,blade_filepath);
                fid0=fopen('cfx_turbogrid.sh','w');
                fprintf(fid0,'%s%d%s%d%s\n','/ansys_inc/v190/TurboGrid/bin/cfxtg -batch /test4T/xiongjin/57fanopt/',56,'/blade',1000*Ge+i+kn-1,'.tse');
                fclose(fid0);
                fprintf('%s%d\n','now running: ',1000*Ge+i+kn-1);
                !sh cfx_turbogrid.sh
                fid4=fopen('cfxpre.sh','w');
                cfxpre(Ge,i+kn-1,N_median,c1(1)/Z(1,1),blade_filepath,blade_filepath)
                fprintf(fid4,'%s%d%s%d%s\n','/ansys_inc/v190/CFX/bin/cfx5pre -batch /test4T/xiongjin/57fanopt/',56,'/blade',1000*Ge+i+kn-1,'.pre');
                fclose(fid4);
                %             cfxpost(Ge,i,blade_filepath)
                pause(5);
                !sh cfxpre.sh
                pause(5);
                !pkill SolverManager
            end
            
            fid1=fopen('cfx_solve.sh','w');
            fprintf(fid1,'%s%d%s\n','echo "',1000*Ge+i+kn-1,' start"');
            fprintf(fid1,'%s%d%s%d%s%d%s\n','/ansys_inc/v190/CFX/bin/cfx5solve -def /test4T/xiongjin/57fanopt/',56,'/blade',1000*Ge+i+kn-1,'.def -part 110 -start-method "Platform MPI Local Parallel" ;');
            fclose(fid1);
            !sh cfx_solve.sh
            AF0=dir(strcat(blade_filepath,'/blade',num2str(1000*Ge+i+kn-1),'*.res'));
            if size(AF0,1)>0
                ABC_final(:,i)=ABC(:,1);
                disp(ABC_final')
                if max(AF0(1).name~=strcat('blade',num2str(1000*Ge+i+kn-1),'_001.res'))~=0
                    movefile(strcat(AF0(1).folder,'/',AF0(1).name),strcat(AF0(1).folder,'/',strcat('blade',num2str(1000*Ge+i+kn-1),'_001.res')));
                    %                     pause(1);
                    %                     delete(strcat(AF0(1).folder,'/',AF0(1).name));
                end
                break;
            end
            clear dCu beta_1A_c beta_1A D_No_blade
        end
        
    end
    
    %=============================
    Af=dir(strcat(blade_filepath,'/blade*.res'));
    Af_done=dir(strcat(blade_filepath,'/output*.txt'));
    
    for i=1:size(Af_done,1)
        res_done(i)=str2num(Af_done(i).name(7:10));
    end
    newcase=0;
    fid1=fopen('cfx_post.sh','w');
    for i=1:size(Af,1)
        res_file(i)=str2num(Af(i).name(6:9));
        if size(find(res_done==res_file(i)),2)==0
            newcase=newcase+1;
            new_res(newcase)=res_file(i);
            cfxpost(res_file(i),blade_filepath);
            fprintf(fid1,'%s%d%s\n','/ansys_inc/v190/CFX/bin/cfx5post -batch /test4T/xiongjin/57fanopt/56/blade',res_file(i),'.cse');
        end
    end
    fclose(fid1);
    !sh cfx_post.sh
    
    if newcase>0
        no_title=1000*Ge+kn-1;
        for i=1:length(new_res)
            file1=sprintf('%s%d%s','output',new_res(i),'.txt');
            Output_data(i,1)=new_res(i);
            Outputdata(i,1:5)=load(file1);
            Output_data(i,2:12)=ABC_final(:,new_res(i)-no_title)';
            c9min=3;c9max=6.999;%<<<<<<<<<<<<<<<<<<c8:cy_est=[0.8 1.8]
            c9=ABC_final(9,new_res(i)-no_title)*(c9max-c9min)+c9min;%
            z=fix(c9);
            Output_data(i,13)=Outputdata(i,5);
            Output_data(i,14)=z;
            Output_data(i,15)=Outputdata(i,1);
            Output_data(i,16)=Outputdata(i,2)*z;
            Output_data(i,17)=Outputdata(i,3);
            Output_data(i,18)=-Outputdata(i,2)/1.185*Outputdata(i,3)/(N_median*Outputdata(i,4)*pi/30)*100;
            Output_data(i,19)=-3600*Outputdata(i,2)/1.185/(N_median*Outputdata(i,4)*pi/30)*100;
            c1min=Qmin;c1max=Qmax;
            c1=ABC_final(1,new_res(i)-no_title)*(c1max-c1min)+c1min;%C1:Qin=[Qmin,Qmax]
            Output_data(i,20)=c1;
            c6min=Pmin;c6max=Pmax;%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            c6=ABC_final(6,new_res(i)-no_title)*(c6max-c6min)+c6min;%C6=[Pmin Pmax]: total pressure rise
            Output_data(i,21)=c6;
            c8min=1.2;c8max=1.8;%<<<<<<<<<<<<<<<<<<c8:cy_est=[0.8 1.8]
            c8=ABC_final(8,new_res(i)-no_title)*(c8max-c8min)+c8min;%
            Output_data(i,22)=c8;
            Output_data(i,23)=Output_data(i,17)/Output_data(i,18)*100;
            Output_data(i,24)=Output_data(i,22)/Output_data(i,21)*Output_data(i,23);
            Output_data(i,25)=-Outputdata(i,2)/1.185*Outputdata(i,5)/(N_median*Outputdata(i,4)*pi/30)*100;
        end
        fid1=fopen('new_result.txt','a');
        for i=1:size(Output_data,1)
            for j=1:size(Output_data,2)
                fprintf(fid1,'%g ',Output_data(i,j));
            end
            fprintf(fid1,'\n');
        end
        fclose(fid1);
        father=[father;[Output_data(:,2:12) Output_data(:,17:18) Output_data(:,13) Output_data(:,25)]];
        flag1=0;
    elseif newcase==0
        if flag1==0
            flag1=1;
        elseif flag1==1
            flag1=2;
        elseif flag1==2
            flag1=3;
        else
            error('converged')
        end
    end
    clear Output_data new_res Outputdata res_file res_done new_offspring new_offspring_b new_offspring_c beta_1A_c beta_1A D_No_blade new_offspring dCu ABC_final
    
end

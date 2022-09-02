function [IGDmean,IGDstd]=compare(fun_name,maxii,pop,gen)
%====================================================================
V=20;%>2
Index=0; %1��ֻ�������ս�� 2������ȫ����� 3������ 4������gif�����ļ�
% ZDT1 ZDT2 ZDT3 ZDT4 ZDT6
% LZ08F1 LZ08F2 LZ08F3 LZ08F4 LZ08F5 LZ08F7 LZ08F8 LZ08F9
% VLMOP2
dx1=0.001;
for ii=1:maxii
%=====================================================================
%% PF of VLMOP2
if strcmp(fun_name,'VLMOP2')
[chromosome,IGD]=nsga_2(pop,gen,2,zeros(1,V),ones(1,V),fun_name,Index,0.001);
if Index~=0
rangex=-1/sqrt(V):dx1:1/sqrt(V);
 for i=1:size(rangex,2)
   for j=1:V
    x(j)=rangex(i);
   end
y0=evaluate_objective(x,fun_name);
range(i)=y0(1);
y(i)=y0(2);   
 end
subplot(1,2,1)
hold on
plot(range,y,'ro')
end
%% PF of ZDT1
elseif strcmp(fun_name,'ZDT1')
[chromosome,IGD]=nsga_2(pop,gen,2,zeros(1,V),ones(1,V),fun_name,Index,0.1);
if Index~=0
range=0:dx1:1;
y=1-sqrt(range);
subplot(1,2,1)
hold on
plot(range,y)
end
elseif strcmp(fun_name,'ZDT2')
%% PF of ZDT2
[chromosome,IGD]=nsga_2(pop,gen,2,zeros(1,V),ones(1,V),fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
y=1-range.^2;
subplot(1,2,1)
hold on
plot(range,y)
end
elseif strcmp(fun_name,'ZDT3')
%% PF of ZDT3
[chromosome,IGD]=nsga_2(pop,gen,2,zeros(1,V),ones(1,V),fun_name,Index,0.001);
if Index~=0
[range,y]=PF(V,dx1,fun_name);
% range=0:dx1:1;
% y=1-sqrt(range)-range.*sin(10*pi*range);
% No=1;
% x0(1)=range(1);
% y0(1)=y(1);
% for i=2:size(range,2)-1
%    if range(i+1)<range(i) || y(i+1)<y(i)
%        No=No+1;
%        x0(No)=range(i+1);
%        y0(No)=y(i);
%    end
% end
subplot(1,2,1)
hold on
plot(range,y,'r*')
end
elseif strcmp(fun_name,'ZDT4')
%% PF of ZDT4
[chromosome,IGD]=nsga_2(pop,gen,2,[0 -5*ones(1,V-1)],[1 5*ones(1,V-1)],fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
y=1-sqrt(range);
subplot(1,2,1)
hold on
plot(range,y)
end
elseif strcmp(fun_name,'ZDT6')
%% PF of ZDT6
[chromosome,IGD]=nsga_2(pop,gen,2,zeros(1,V),ones(1,V),fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
y=1-range.^2;
subplot(1,2,1)
hold on
plot(range,y)
end
elseif strcmp(fun_name,'LZ08F1')
%% PF of LZ08-F1
[chromosome,IGD]=nsga_2(pop,gen,2,zeros(1,V),ones(1,V),fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
for j=1:size(range,2)
    x(1)=range(j);
for i=2:20
    x(i)=x(1)^(0.5+1.5*(i-2)/18);
end
y=evaluate_objective(x,fun_name);
subplot(1,2,1)
hold on
plot(y(1),y(2),'ro')
x1(j)=x(1);
x2(j)=x(2);
x3(j)=x(3);
end
figure(2)
plot3(chromosome(:,1)',chromosome(:,2)',chromosome(:,3)','ro')
hold on
plot3(x1,x2,x3,'g*')
end
elseif strcmp(fun_name,'LZ08F2')
%% PF of LZ08-F2
[chromosome,IGD]=nsga_2(pop,gen,2,[0 -ones(1,V-1)],[1 ones(1,V-1)],fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
for j=1:size(range,2)
    x(1)=range(j);
for i=2:V
    x(i)=sin(6*pi*x(1)+i*pi/V);
end
y=evaluate_objective(x,fun_name);
subplot(1,2,1)
hold on
plot(y(1),y(2),'ro')
x1(j)=x(1);
x2(j)=x(2);
x3(j)=x(3);
end
figure(2)
plot3(chromosome(:,1)',chromosome(:,2)',chromosome(:,3)','ro')
hold on
plot3(x1,x2,x3,'g*')
end
elseif strcmp(fun_name,'LZ08F3')
%% PF of LZ08-F3
[chromosome,IGD]=nsga_2(pop,gen,2,[0 -ones(1,V-1)],[1 ones(1,V-1)],fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
for j=1:size(range,2)
    x(1)=range(j);
    
for i=2:V
    if floor(i/2)==i/2
    x(i)=0.8*x(1)*sin(6*pi*x(1)+i*pi/V);
    else
    x(i)=0.8*x(1)*cos(6*pi*x(1)+i*pi/V);
    end
end

y=evaluate_objective(x,fun_name);
subplot(1,2,1)
hold on
plot(y(1),y(2),'ro')
x1(j)=x(1);
x2(j)=x(2);
x3(j)=x(3);
end
figure(2)
plot3(chromosome(:,1)',chromosome(:,2)',chromosome(:,3)','ro')
hold on
plot3(x1,x2,x3,'g*')
end
elseif strcmp(fun_name,'LZ08F4')
%% PF of LZ08-F4
[chromosome,IGD]=nsga_2(pop,gen,2,[0 -ones(1,V-1)],[1 ones(1,V-1)],fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
for j=1:size(range,2)
    x(1)=range(j);
    
for i=2:V
    if floor(i/2)==i/2
    x(i)=0.8*x(1)*sin(6*pi*x(1)+i*pi/V);
    else
    x(i)=0.8*x(1)*cos((6*pi*x(1)+i*pi/V)/3);
    end
end

y=evaluate_objective(x,fun_name);
subplot(1,2,1)
hold on
plot(y(1),y(2),'ro')
x1(j)=x(1);
x2(j)=x(2);
x3(j)=x(3);
end
figure(2)
plot3(chromosome(:,1)',chromosome(:,2)',chromosome(:,3)','ro')
hold on
plot3(x1,x2,x3,'g*')
end
elseif strcmp(fun_name,'LZ08F5')
%% PF of LZ08-F5
[chromosome,IGD]=nsga_2(pop,gen,2,[0 -ones(1,V-1)],[1 ones(1,V-1)],fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
for j=1:size(range,2)
    x(1)=range(j);
    
for i=2:V
    if floor(i/2)==i/2
    x(i)=(0.3*x(1)^2*cos(24*pi*x(1)+4*i*pi/V)+0.6*x(1))*sin(6*pi*x(1)+i*pi/V);
    else
    x(i)=(0.3*x(1)^2*cos(24*pi*x(1)+4*i*pi/V)+0.6*x(1))*cos(6*pi*x(1)+i*pi/V);
    end
end

y=evaluate_objective(x,fun_name);
subplot(1,2,1)
hold on
plot(y(1),y(2),'ro')
x1(j)=x(1);
x2(j)=x(2);
x3(j)=x(3);
end
figure(2)
plot3(chromosome(:,1)',chromosome(:,2)',chromosome(:,3)','ro')
hold on
plot3(x1,x2,x3,'g*')
end
elseif strcmp(fun_name,'LZ08F7')||strcmp(fun_name,'LZ08F8')
%% PF of LZ08-F7 and LZ08-F8
[chromosome,IGD]=nsga_2(pop,gen,2,zeros(1,V),ones(1,V),fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
for j=1:size(range,2)
    x(1)=range(j);
    
for i=2:V
    x(i)=x(1)^(0.5*(1+3*(i-2)/(V-2)));
end

y=evaluate_objective(x,fun_name);
subplot(1,2,1)
hold on
plot(y(1),y(2),'ro')
x1(j)=x(1);
x2(j)=x(2);
x3(j)=x(3);
end
figure(2)
plot3(chromosome(:,1)',chromosome(:,2)',chromosome(:,3)','ro')
hold on
plot3(x1,x2,x3,'g*')
end
elseif strcmp(fun_name,'LZ08F9')
%% PF of LZ08-F9
[chromosome,IGD]=nsga_2(pop,gen,2,[0 -ones(1,V-1)],[1 ones(1,V-1)],fun_name,Index,0.001);
if Index~=0
range=0:dx1:1;
for j=1:size(range,2)
    x(1)=range(j); 
for i=2:V
    x(i)=sin(6*pi*x(1)+i*pi/V);
end
y=evaluate_objective(x,fun_name);
subplot(1,2,1)
hold on
plot(y(1),y(2),'ro')
x1(j)=x(1);
x2(j)=x(2);
x3(j)=x(3);
end
figure(2)
plot3(chromosome(:,1)',chromosome(:,2)',chromosome(:,3)','ro')
hold on
plot3(x1,x2,x3,'g*')
end
else
error('�����������еĺ���')     
end
IGDii(ii)=IGD(gen);
end
IGDmean=mean(IGDii);
IGDstd=std(IGDii);
fileID = fopen('log.txt','a+');
fprintf(fileID , '%s%d%s%d%s%s%s%d%s%2.6f%s%2.6f\n' , 'pop=',pop,' gen=',gen,' function:',fun_name,' ii=',maxii,' IGD=',IGDmean,'��',IGDstd);
fclose(fileID);
end
function [range,y]=PF(V,dx1,fun_name)
if strcmp(fun_name,'VLMOP2')
%% PF of VLMOP2
rangex=-1/sqrt(V):dx1/2:1/sqrt(V);
for i=1:size(rangex,2)
   for j=1:V
    x(j)=rangex(i);
   end
y0=evaluate_objective(x,fun_name);
range(i)=y0(1);
y(i)=y0(2);   
end
elseif strcmp(fun_name,'ZDT1')
%% PF of ZDT1
range=[0:dx1:1];
y=1-sqrt(range);
elseif strcmp(fun_name,'ZDT2')
%% PF of ZDT2
range=[0:dx1:1];
y=1-range.^2;
elseif strcmp(fun_name,'ZDT3')
%% PF of ZDT3
range0=[0:dx1/4:1];
y0=1-sqrt(range0)-range0.*sin(10*pi*range0);
Jud=1;
No=0;
for i=1:size(range0,2)
    for j=1:size(range0,2)
      if (range0(i)>=range0(j) && y0(i)>=y0(j) ) && ~(range0(i)==range0(j) && y0(i)==y0(j))  
       Jud=0;
       break
      end
    end
if Jud
  No=No+1;
  range(No)=range0(i);
  y(No)=y0(i);
end
Jud=1;    
end
elseif strcmp(fun_name,'ZDT4')
%% PF of ZDT4 
range=[0:dx1:1];
y=1-sqrt(range);
elseif strcmp(fun_name,'ZDT6')
%% PF of ZDT6
range=[0.281:dx1*0.72:1];
y=1-range.^2;
elseif strcmp(fun_name,'LZ08F1')
%% PF of LZ08-F1
range=[0:dx1:1];
for j=1:size(range,2)
    x(1)=range(j);
for i=2:20
    x(i)=x(1)^(0.5+1.5*(i-2)/18);
end
y0=evaluate_objective(x,fun_name);
range(j)=y0(1);
y(j)=y0(2);
end
elseif strcmp(fun_name,'LZ08F2')
%% PF of LZ08-F2

range=[0:dx1:1];
for j=1:size(range,2)
    x(1)=range(j);
for i=2:V
    x(i)=sin(6*pi*x(1)+i*pi/V);
end
y0=evaluate_objective(x,fun_name);
range(j)=y0(1);
y(j)=y0(2);
end
elseif strcmp(fun_name,'LZ08F3')
%% PF of LZ08-F3

range=[0:dx1:1];
for j=1:size(range,2)
    x(1)=range(j);
    
for i=2:V
    if floor(i/2)==i/2
    x(i)=0.8*x(1)*sin(6*pi*x(1)+i*pi/V);
    else
    x(i)=0.8*x(1)*cos(6*pi*x(1)+i*pi/V);
    end
end
y0=evaluate_objective(x,fun_name);
range(j)=y0(1);
y(j)=y0(2);
end
elseif strcmp(fun_name,'LZ08F4')
%% PF of LZ08-F4
range=[0:dx1:1];
for j=1:size(range,2)
    x(1)=range(j);
    
for i=2:V
    if floor(i/2)==i/2
    x(i)=0.8*x(1)*sin(6*pi*x(1)+i*pi/V);
    else
    x(i)=0.8*x(1)*cos((6*pi*x(1)+i*pi/V)/3);
    end
end
y0=evaluate_objective(x,fun_name);
range(j)=y0(1);
y(j)=y0(2);
end
elseif strcmp(fun_name,'LZ08F5')
%% PF of LZ08-F5

range=[0:dx1:1];
for j=1:size(range,2)
    x(1)=range(j);
    
for i=2:V
    if floor(i/2)==i/2
    x(i)=(0.3*x(1)^2*cos(24*pi*x(1)+4*i*pi/V)+0.6*x(1))*sin(6*pi*x(1)+i*pi/V);
    else
    x(i)=(0.3*x(1)^2*cos(24*pi*x(1)+4*i*pi/V)+0.6*x(1))*cos(6*pi*x(1)+i*pi/V);
    end
end
y0=evaluate_objective(x,fun_name);
range(j)=y0(1);
y(j)=y0(2);
end
elseif strcmp(fun_name,'LZ08F7')||strcmp(fun_name,'LZ08F8')
%% PF of LZ08-F7 and LZ08-F8

range=[0:dx1:1];
for j=1:size(range,2)
    x(1)=range(j);
    
for i=2:V
    x(i)=x(1)^(0.5*(1+3*(i-2)/(V-2)));
end
y0=evaluate_objective(x,fun_name);
range(j)=y0(1);
y(j)=y0(2);
end
elseif strcmp(fun_name,'LZ08F9')
%% PF of LZ08-F9
range=[0:dx1:1];
for j=1:size(range,2)
    x(1)=range(j); 
for i=2:V
    x(i)=sin(6*pi*x(1)+i*pi/V);
end
y0=evaluate_objective(x,fun_name);
range(j)=y0(1);
y(j)=y0(2);
end
else
error('仅限输入已有的函数')     
end

function f = evaluate_objective(x,name)
f = [];
V=size(x,2);
J1=floor(V/2-0.5);
J2=floor(V/2);
if strcmp(name,'ZDT1')
%% Deb, Kalyanmoy. "Multi-objective genetic algorithms: Problem difficulties and construction of test problems." Evolutionary computation 7.3 (1999): 205-230.
% % ZDT1
% % y1,y2,xi=[0,1]^N 
f(1) = x(1);
g = 1;
sum0 = 0;
for i = 2:V
    sum0 = sum0 + x(i);
end
sum0 = 9*(sum0 / (V-1));
g = g + sum0;
f(2) = g * (1 - sqrt(x(1) / g));
elseif strcmp(name,'ZDT2')
%% Deb, Kalyanmoy. "Multi-objective genetic algorithms: Problem difficulties and construction of test problems." Evolutionary computation 7.3 (1999): 205-230.
% % ZDT2
% % y1,y2,xi=[0,1]^N
f(1) = x(1);
g = 1;
sum0 = 0;
for i = 2:V
    sum0 = sum0 + x(i);
end
sum0 = 9*(sum0 / (V-1));
g = g + sum0;
f(2) = g * (1 - (x(1) / g)^2);
elseif strcmp(name,'ZDT3')
%% Deb, Kalyanmoy. "Multi-objective genetic algorithms: Problem difficulties and construction of test problems." Evolutionary computation 7.3 (1999): 205-230.
% % ZDT3
% % y1,y2,xi=[0,1]^N 
f(1) = x(1);
g = 1;
sum0 = 0;
for i = 2:V
    sum0 = sum0 + x(i);
end
sum0 = 9*(sum0 / (V-1));
g = g + sum0;
f(2) = g * (1 - sqrt(x(1) / g)-x(1)/g*sin(10*pi*x(1)));
elseif strcmp(name,'ZDT4')
%% Deb, Kalyanmoy. "Multi-objective genetic algorithms: Problem difficulties and construction of test problems." Evolutionary computation 7.3 (1999): 205-230.
% % ZDT4
% % y1,y2,x=[0,1]x[-5,5]^N-1
f(1) = x(1);
g = 1+10*(V-1);
for i = 2:V
    g = g + x(i)^2-10*cos(4*pi*x(i));
end
f(2) = g * (1 - sqrt(x(1) / g));
elseif strcmp(name,'ZDT6')
%% Deb, Kalyanmoy. "Multi-objective genetic algorithms: Problem difficulties and construction of test problems." Evolutionary computation 7.3 (1999): 205-230.
% % ZDT6
% % y1,y2,xi=[0,1]^N 
f(1) = 1-exp(-4*x(1))*(sin(6*pi*x(1)))^6;
g = 1;
sum0=0;
for i = 2:V
    sum0 = sum0 + x(i);
end
g = g + 9*(sum0/(V-1))^0.25;
f(2) = g * (1 - ( f(1) / g)^2);
elseif strcmp(name,'VLMOP2')
%% Van Veldhuizen, David A., and Gary B. Lamont. "Multiobjective evolutionary algorithm test suites." SAC. Vol. 99. 1999.
% VLMOP2
% y1,y2,x=[-2,2]^N
sum1 = 0;
sum2 = 0;
for i=1:V
  sum1=sum1+(x(i)-1/sqrt(V))^2;
  sum2=sum2+(x(i)+1/sqrt(V))^2;
end
f(1) =1-exp(-sum1) ;
f(2) =1-exp(-sum2) ;
elseif strcmp(name,'LZ08F1')
%% Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." IEEE transactions on evolutionary computation 13.2 (2008): 284-302.
% LZ08-F1
% y1,y2,x=[0,1]^N
odd=0;
even=0;
for i=2:V
    if floor(i/2)==i/2
        even=even+2/J2*(x(i)-x(1)^(0.5+1.5*(i-2)/(V-2)))^2;
    else
        odd=odd+2/J1*(x(i)-x(1)^(0.5+1.5*(i-2)/(V-2)))^2;
    end
end        
f(1) = x(1)+odd;
f(2) = 1-sqrt(x(1))+even;
elseif strcmp(name,'LZ08F2')
%% Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." IEEE transactions on evolutionary computation 13.2 (2008): 284-302.
% LZ08-F2
% y1,y2,x=[0,1]x[-1,1]^N-1
odd=0;
even=0;
for i=2:V
    if floor(i/2)==i/2
        even=even+2/J2*(x(i)-sin(6*pi*x(1)+i*pi/V))^2;
    else
        odd = odd+2/J1*(x(i)-sin(6*pi*x(1)+i*pi/V))^2;
    end
end 
f(1) = x(1)+odd;
f(2) = 1-sqrt(x(1))+even;
elseif strcmp(name,'LZ08F3')
%% Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." IEEE transactions on evolutionary computation 13.2 (2008): 284-302.
% LZ08-F3
% y1,y2,x=[0,1]x[-1,1]^N-1
odd=0;
even=0;
for i=2:V
    if floor(i/2)==i/2
        even=even+2/J2*(x(i)-0.8*x(1)*sin(6*pi*x(1)+i*pi/V))^2;
    else
        odd = odd+2/J1*(x(i)-0.8*x(1)*cos(6*pi*x(1)+i*pi/V))^2;
    end
end 
f(1) = x(1)+odd;
f(2) = 1-sqrt(x(1))+even;
elseif strcmp(name,'LZ08F4')
%% Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." IEEE transactions on evolutionary computation 13.2 (2008): 284-302.
% LZ08-F4
% y1,y2,x=[0,1]x[-1,1]^N-1
odd=0;
even=0;
for i=2:V
    if floor(i/2)==i/2
        even=even+2/J2*(x(i)-0.8*x(1)*sin(6*pi*x(1)+i*pi/V))^2;
    else
        odd = odd+2/J1*(x(i)-0.8*x(1)*cos((6*pi*x(1)+i*pi/V)/3))^2;
    end
end 
f(1) = x(1)+odd;
f(2) = 1-sqrt(x(1))+even;
elseif strcmp(name,'LZ08F5')
%% Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." IEEE transactions on evolutionary computation 13.2 (2008): 284-302.
% LZ08-F5
% y1,y2,x=[0,1]x[-1,1]^N-1
odd=0;
even=0;
for i=2:V
    if floor(i/2)==i/2
        even=even+2/J2*(x(i)-(0.3*x(1)^2*cos(24*pi*x(1)+4*i*pi/V)+0.6*x(1))*sin(6*pi*x(1)+i*pi/V))^2;
    else
        odd = odd+2/J1*(x(i)-(0.3*x(1)^2*cos(24*pi*x(1)+4*i*pi/V)+0.6*x(1))*cos(6*pi*x(1)+i*pi/V))^2;
    end
end 
f(1) = x(1)+odd;
f(2) = 1-sqrt(x(1))+even;
elseif strcmp(name,'LZ08F7')
%% Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." IEEE transactions on evolutionary computation 13.2 (2008): 284-302.
% LZ08-F7
% y1,y2,x=[0,1]^N
odd=0;
even=0;
for i=2:V
    y0=x(i)-x(1)^(0.5*(1+3*(i-2)/(V-2)));
    if floor(i/2)==i/2
        even=even+2/J2*(4*y0^2-cos(8*y0*pi)+1);
    else
        odd = odd+2/J1*(4*y0^2-cos(8*y0*pi)+1);
    end
end 
f(1) = x(1)+odd;
f(2) = 1-sqrt(x(1))+even;
elseif strcmp(name,'LZ08F8')
%% Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." IEEE transactions on evolutionary computation 13.2 (2008): 284-302.
% LZ08-F8
% y1,y2,x=[0,1]^N
odd1=0;
odd2=1;
even1=0;
even2=1;
for i=2:V
    y0=x(i)-x(1)^(0.5*(1+3*(i-2)/(V-2)));
    if floor(i/2)==i/2
        even1=even1+y0^2;
        even2=even2*cos(20*y0*pi/sqrt(i));
    else
        odd1 = odd1+y0^2;
        odd2 = odd2*cos(20*y0*pi/sqrt(i));
    end
end 
f(1) = x(1)+2/J1*(4*odd1-2*odd2+2);
f(2) = 1-sqrt(x(1))+2/J2*(4*odd1-2*odd2+2);
elseif strcmp(name,'LZ08F9')
%% Li, Hui, and Qingfu Zhang. "Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II." IEEE transactions on evolutionary computation 13.2 (2008): 284-302.
% LZ08-F9
% y1,y2,x=[0,1]x[-1,1]^N-1
odd=0;
even=0;
for i=2:V
    if floor(i/2)==i/2
        even=even+2/J2*(x(i)-sin(6*pi*x(1)+i*pi/V))^2;
    else
        odd = odd+2/J1*(x(i)-sin(6*pi*x(1)+i*pi/V))^2;
    end
end 
f(1) = x(1)+odd;
f(2) = 1-x(1)^2+even;
else
error('仅限输入已有的函数') 
end
end

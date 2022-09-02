function meanline_xy=meanline_parabola(b,beta_1A,beta_2A,ratio,number)
if beta_2A>beta_1A
    beta_m=(beta_2A-beta_1A)*(1-ratio)+beta_1A;
    x1=beta_m-beta_1A;
    x2=beta_2A-beta_m;
    A=(cotd(x2)-cotd(x1))/2;
    B=-b;
    C=b*cotd(x1);
    x=linspace(0,b,number);
    y(1)=0;y(number)=0;
    if A~=0
        for i=2:length(x)-1
            a0=A^2;
            b0=2*A*x(i)+C;
            c0=x(i)^2+B*x(i);
            %     y(i)=max(solve(@(y0)x(i)^2+2*A*x(i)*y0+A^2*y0^2+B*x(i)+C*y0==0));
            y(i)=(-b0+(b0^2-4*a0*c0)^0.5)/(2*a0);
        end  
    elseif A==0
        for i=2:length(x)-1
            y(i)=-(x(i)^2+B*x(i))/C;
        end
    end
    x=linspace(-b/2,b/2,number);
    x0=x*cosd(beta_m)+y*sind(beta_m);
    y0=y*cosd(beta_m)-x*sind(beta_m);
%     plot(x,y,'r')
%     hold on
%     plot(x0,y0,'b')
    meanline_xy=[x0;y0];
    
elseif beta_2A<beta_1A
    beta_m=(beta_2A-beta_1A)*ratio+beta_1A;
    x2=-(beta_m-beta_1A);
    x1=-(beta_2A-beta_m);
    A=(cotd(x2)-cotd(x1))/2;
    B=-b;
    C=b*cotd(x1);
    x=linspace(0,b,number);
    y(1)=0;y(number)=0;
    if A~=0
        for i=2:length(x)-1
            a0=A^2;
            b0=2*A*x(i)+C;
            c0=x(i)^2+B*x(i);
            %     y(i)=max(solve(@(y0)x(i)^2+2*A*x(i)*y0+A^2*y0^2+B*x(i)+C*y0==0));
            y(i)=-(-b0+(b0^2-4*a0*c0)^0.5)/(2*a0);
        end  
    elseif A==0
        for i=2:length(x)-1
            y(i)=-(x(i)^2+B*x(i))/C;
        end
    end
    x=linspace(-b/2,b/2,number);
    x0=x*cosd(beta_1A-x1)+y*sind(beta_1A-x1);
    y0=y*cosd(beta_1A-x1)-x*sind(beta_1A-x1);
%     plot(x,y,'r')
%     hold on
%     plot(x0,y0,'b')
    meanline_xy=[x0;y0];   
else
    x=linspace(-b/2*cosd(beta_2A),b/2*cosd(beta_2A),number);
    y=linspace(-b/2*sind(-beta_2A),b/2*sind(-beta_2A),number);
%     plot(x,y,'r')
    meanline_xy=[x;y];
end
end
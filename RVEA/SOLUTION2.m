function [obj,FE] = SOLUTION2(PopDec,FE,ProblemM,index)
%可以换成任意目标函数
FE=FE+size(PopDec,1);
ProblemD=size(PopDec,2);
switch index
    case 'test'
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            for j=1:ProblemM
                obj(i).obj(j)=sum(PopDec(i,1:end-2*(j-1)))-sum(PopDec(i,end-2*(j-1):end));
            end
        end

    case 'DTLZ1'
        g   =100*(ProblemD-ProblemM+1+sum((PopDec(:,ProblemM:end)-0.5).^2-cos(20.*pi.*(PopDec(:,ProblemM:end)-0.5)),2));
        objobj = 0.5*repmat(1+g,1,ProblemM).*fliplr(cumprod([ones(size(PopDec,1),1),PopDec(:,1:ProblemM-1)],2)).*[ones(size(PopDec,1),1),1-PopDec(:,ProblemM-1:-1:1)];
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'DTLZ2'
        g      = sum((PopDec(:,ProblemM:end)-0.5).^2,2);
        objobj = repmat(1+g,1,ProblemM).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:ProblemM-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,ProblemM-1:-1:1)*pi/2)];
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'DTLZ3'
        g      = 100*(ProblemD-ProblemM+1+sum((PopDec(:,ProblemM:end)-0.5).^2-cos(20.*pi.*(PopDec(:,ProblemM:end)-0.5)),2));
        objobj = repmat(1+g,1,ProblemM).*fliplr(cumprod([ones(size(PopDec,1),1),cos(PopDec(:,1:ProblemM-1)*pi/2)],2)).*[ones(size(PopDec,1),1),sin(PopDec(:,ProblemM-1:-1:1)*pi/2)];
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'DTLZ4'
        PopDec(:,1:ProblemM-1) = PopDec(:,1:ProblemM-1).^100;
        g      = sum((PopDec(:,ProblemM:end)-0.5).^2,2);
        objobj = repmat(1+g,1,ProblemM).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:ProblemM-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,ProblemM-1:-1:1)*pi/2)];
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'DTLZ5'
        g    = sum((PopDec(:,ProblemM:end)-0.5).^2,2);
        Temp = repmat(g,1,ProblemM-2);
        PopDec(:,2:ProblemM-1) = (1+2*Temp.*PopDec(:,2:ProblemM-1))./(2+2*Temp);
        objobj = repmat(1+g,1,ProblemM).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:ProblemM-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,ProblemM-1:-1:1)*pi/2)];
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'DTLZ6'
        g    = sum(PopDec(:,ProblemM:end).^0.1,2);
        Temp = repmat(g,1,ProblemM-2);
        PopDec(:,2:ProblemM-1) = (1+2*Temp.*PopDec(:,2:ProblemM-1))./(2+2*Temp);
        objobj = repmat(1+g,1,ProblemM).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:ProblemM-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,ProblemM-1:-1:1)*pi/2)];
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end


    case 'DTLZ7'
        objobj = zeros(size(PopDec,1),ProblemM);
        g      = 1+9*mean(PopDec(:,ProblemM:end),2);
        objobj(:,1:ProblemM-1) = PopDec(:,1:ProblemM-1);
        objobj(:,ProblemM)     = (1+g).*(ProblemM-sum(objobj(:,1:ProblemM-1)./(1+repmat(g,1,ProblemM-1)).*(1+sin(3*pi.*objobj(:,1:ProblemM-1))),2));
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'ZDT1'
        if ProblemM~=2
            error('error');
        end
        objobj(:,1) = PopDec(:,1);
        g = 1 + 9*mean(PopDec(:,2:end),2);
        h = 1 - (objobj(:,1)./g).^0.5;
        objobj(:,2) = g.*h;
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'ZDT2'
        if ProblemM~=2
            error('error');
        end
        objobj(:,1) = PopDec(:,1);
        g = 1 + 9*mean(PopDec(:,2:end),2);
        h = 1 - (objobj(:,1)./g).^2;
        objobj(:,2) = g.*h;
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'ZDT3'
        if ProblemM~=2
            error('error');
        end
        objobj(:,1) = PopDec(:,1);
        g = 1 + 9*mean(PopDec(:,2:end),2);
        h = 1 - (objobj(:,1)./g).^0.5 - objobj(:,1)./g.*sin(10*pi*objobj(:,1));
        objobj(:,2) = g.*h;
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'ZDT4'
        if ProblemM~=2
            error('error');
        end
        objobj(:,1) = PopDec(:,1);
        g = 1 + 10*(size(PopDec,2)-1) + sum(PopDec(:,2:end).^2-10*cos(4*pi*PopDec(:,2:end)),2);
        h = 1 - (objobj(:,1)./g).^0.5;
        objobj(:,2) = g.*h;
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'ZDT5'
        if ProblemM~=2
            error('error');
        end
        u      = zeros(size(PopDec,1),1+(size(PopDec,2)-30)/5);
        u(:,1) = sum(PopDec(:,1:30),2);
        for i = 2 : size(u,2)
            u(:,i) = sum(PopDec(:,(i-2)*5+31:(i-2)*5+35),2);
        end
        v           = zeros(size(u));
        v(u<5)      = 2 + u(u<5);
        v(u==5)     = 1;
        objobj(:,1) = 1 + u(:,1);
        g           = sum(v(:,2:end),2);
        h           = 1./objobj(:,1);
        objobj(:,2) = g.*h;
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

    case 'ZDT6'
        if ProblemM~=2
            error('error');
        end
        objobj(:,1) = 1 - exp(-4*PopDec(:,1)).*sin(6*pi*PopDec(:,1)).^6;
        g = 1 + 9*mean(PopDec(:,2:end),2).^0.25;
        h = 1 - (objobj(:,1)./g).^2;
        objobj(:,2) = g.*h;
        for i=1:size(PopDec,1)
            obj(i).dec=PopDec(i,:);
            obj(i).obj=objobj(i,:);
        end

end
end
function R=PF(obj,N)

switch obj.name
    case 'DTLZ1'
        R = UniformPoint2(N,obj.M)/2;
        %         if Problem.M==2
        %             R=UniformPoint2(N,Problem.M)/2;
        %         elseif Problem.M==3
        %             a = linspace(0,1,10)';
        %             R = {a*a'/2,a*(1-a')/2,(1-a)*ones(size(a'))/2};
        %         else
        %             R=[];
        %         end
    case 'DTLZ2'
        R = UniformPoint2(N,obj.M);
        R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);

    case 'DTLZ3'
        R = UniformPoint2(N,obj.M);
        R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);

    case 'DTLZ4'
        R = UniformPoint2(N,obj.M);
        R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);

    case 'DTLZ5'
        R = [0:1/(N-1):1;1:-1/(N-1):0]';
        R = R./repmat(sqrt(sum(R.^2,2)),1,size(R,2));
        R = [R(:,ones(1,obj.M-2)),R];
        R = R./sqrt(2).^repmat([obj.M-2,obj.M-2:-1:0],size(R,1),1);

    case 'DTLZ6'
        R = [0:1/(N-1):1;1:-1/(N-1):0]';
        R = R./repmat(sqrt(sum(R.^2,2)),1,size(R,2));
        R = [R(:,ones(1,obj.M-2)),R];
        R = R./sqrt(2).^repmat([obj.M-2,obj.M-2:-1:0],size(R,1),1);

    case 'DTLZ7'
        interval     = [0,0.251412,0.631627,0.859401];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
        X            = UniformPoint2(N,obj.M-1,'grid');
        X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        R            = [X,2*(obj.M-sum(X/2.*(1+sin(3*pi.*X)),2))];

    case 'ZDT1'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^0.5;

    case 'ZDT2'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;

    case 'ZDT3'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^0.5 - R(:,1).*sin(10*pi*R(:,1));
        R      = R(NDSort(R,1)==1,:);

    case 'ZDT4'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^0.5;

    case 'ZDT5'
        R(:,1) = 1 : 31;
        R(:,2) = (obj.D-30)./5./R(:,1);
        
    case 'ZDT6'
        minf1  = 0.280775;
        R(:,1) = linspace(minf1,1,N)';
        R(:,2) = 1 - R(:,1).^2;
end
end
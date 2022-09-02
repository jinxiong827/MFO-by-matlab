function [chromosome,IGD]=nsga_2(pop,gen, M, min_range,max_range,name,Index,dt)

% pop - Population size
% gen - Total number of generations
% M - number of objective function 
% min_range - minimum of the variables e.g.[0 0 0 0]
% max_range - maximum of the variables e.g.[1 1 1 1]
% name - function name in evaluate_objective.m
% Make sure pop and gen are integers
pop = round(pop);
gen = round(gen);
Vmin = size (min_range,2);
Vmax = size (max_range,2);
if Vmin==Vmax
   V=Vmin;
else
    error('��ȷ�����б������ϡ�����')
end
%% Initialize the population
% Population is initialized with random values which are within the
% specified range. Each chromosome consists of the decision variables. Also
% the value of the objective functions, rank and crowding distance
% information is also added to the chromosome vector but only the elements
% of the vector which has the decision variables are operated upon to
% perform the genetic operations like corssover and mutation.
chromosome = initialize_variables(pop, M, V, min_range, max_range,name);


%% Sort the initialized population
% Sort the population using non-domination-sort. This returns two columns
% for each individual which are the rank and the crowding distance
% corresponding to their position in the front they belong. At this stage
% the rank and the crowding distance for each chromosome is added to the
% chromosome vector for easy of computation.
chromosome = non_domination_sort_mod(chromosome, M, V);

%% Start the evolution process
% The following are performed in each generation
% * Select the parents which are fit for reproduction
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.
pic_num = 1;

pn_train=[];
tn_train1=[];
tn_train2=[];
theta = 5*ones(1,20); lob = 0.1*ones(1,20); upb = 20*ones(1,20);
for i = 1 : gen
    % Select the parents
    % Parents are selected for reproduction to generate offspring. The
    % original NSGA-II uses a binary tournament selection based on the
    % crowded-comparision operator. The arguments are 
    % pool - size of the mating pool. It is common to have this to be half the
    %        population size.
    % tour - Tournament size. Original NSGA-II uses a binary tournament
    %        selection, but to see the effect of tournament size this is kept
    %        arbitary, to be choosen by the user.
    pool = round(pop/2);
    tour = 2;
    % Selection process
%     parent_chromosome = tournament_selection(chromosome, pool, tour);
    ind_test=0;
    ind_train=size(pn_train,1);
    pn_test=[];tn_test1=[];tn_test2=[];
    for j=1:pop
        if chromosome(j,23)==1
            ind_test=ind_test+1;
            parent_chromosome0(ind_test,:)=chromosome(j,:);
            pn_test(ind_test,:)=chromosome(j,1:V);
            tn_test1(ind_test)=chromosome(j,V+1);
            tn_test2(ind_test)=chromosome(j,V+2);
        else
            ind_train=ind_train+1;
            pn_train(ind_train,:)=chromosome(j,1:V);
            tn_train1(ind_train)=chromosome(j,V+1);
            tn_train2(ind_train)=chromosome(j,V+2);                     
        end
    end

    %ѵ��SVM
    [mS,mY]=dsmerge(chromosome(:,1:V),chromosome(:,V+1:V+2));
    [dmodel1, perf] =dacefit(mS, mY(:,1), @regpoly1, @corrgauss, theta, lob, upb);
    [dmodel2, perf] =dacefit(mS, mY(:,2), @regpoly1, @corrgauss, theta, lob, upb);      
    save modeldata1 dmodel1 
    save modeldata2 dmodel2 
       x00=[];
    fval=[];
    fitnessfcn=@object;
    %�ҵ�SVMģ�͵�paretoǰ�ؽ�
    options=gaoptimset('paretoFraction',pool/100,'populationsize',100,'generations',...
    500,'stallGenLimit',200,'TolFun',1e-100);
%     options=gaoptimset('paretoFraction',pool/100,'populationsize',100,'generations',...
%     500,'stallGenLimit',200,'TolFun',1e-100,'PlotFcns',@gaplotpareto);
    [x00,fval]=gamultiobj(fitnessfcn,V,[],[],[],[],min_range,max_range,options);
    
%=============================================================================    
    guide_parent = sortrows(non_domination_sort_mod([x00 fval], M, V),V+M+2);
    parent_chromosome=sortrows(parent_chromosome0,V+M+2);
    pn_train((i-1)*pop+1:i*pop,:)=chromosome(:,1:V);
    tn_train1((i-1)*pop+1:i*pop)=chromosome(:,V+1);
    tn_train2((i-1)*pop+1:i*pop)=chromosome(:,V+2); 
    
    [aa1,~]=size(guide_parent);
    [aa2,~]=size(parent_chromosome);
    guide_judgement0=[];
    guide_judgement=[];
    
    for m=1:aa1
        for n=1:V+2
            if n==1
            guide_judgement0(m,n)=0;
            else
            guide_judgement0(m,n)=guide_parent(m,n);
            end
        end
    end
    
    for m=aa1+1:aa1+aa2
        for n=1:V+2
            if n==1
            guide_judgement0(m,n)=1;
            else
            guide_judgement0(m,n)=parent_chromosome(m-aa1,n);
            end 
        end
    end 
    
    guide_judgement=non_domination_sort_mod(guide_judgement0, M, V);
    judgement=guide_judgement(1:aa1,1);
%=============================================================================      
%%  
    mu = 20;
    mum = 20;
    if sum(judgement)<10
    offspring_chromosome = ...
        genetic_operator(parent_chromosome,  guide_parent,...
        pop, V, mu, mum, min_range, max_range);
    else
    offspring_chromosome = ...
        genetic_operator_original(parent_chromosome,...
        pop, V, mu, mum, min_range, max_range);
    end
     judgementrecord(i)=sum(judgement);

    % Intermediate population
    % Intermediate population is the combined population of parents and
    % offsprings of the current generation. The population size is two
    % times the initial population.
    
    [main_pop,~] = size(chromosome);
    [offspring_pop,~] = size(offspring_chromosome);
    % intermediate_chromosome is a concatenation of current population and
    % the offspring population.
    maxi=size(offspring_chromosome,1);
    for j = 1 : maxi
      offspring_chromosome(j,V + 1: M + V) = evaluate_objective(offspring_chromosome(j,:),name);
    end
    intermediate_chromosome(1:main_pop,:) = chromosome;  
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = ...
        offspring_chromosome;

    % Non-domination-sort of intermediate population
    % The intermediate population is sorted again based on non-domination sort
    % before the replacement operator is performed on the intermediate
    % population.
    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, M, V);
    % Perform Selection
    % Once the intermediate population is sorted only the best solution is
    % selected based on it rank and crowding distance. Each front is filled in
    % ascending order until the addition of population size is reached. The
    % last front is included in the population based on the individuals with
    % least crowding distance
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
    fileID0 = fopen('process.txt','w');
    fprintf(fileID0 , '%s\n%d',name,i);
    fclose(fileID0);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%
    [range,y]=PF(V,0.001,name);
    No1=size(find(chromosome(:,23)==1),1);
    PFno1=chromosome(1:No1,21:22);
    distance0=zeros(1,No1);
    for j=1:No1
        distance0(j)=min((PFno1(j,1)-range).^2+(PFno1(j,2)-y).^2);       
    end
    
    IGD(i)=mean(distance0.^0.5);
%     if ~mod(i,10)
%         clc
%         fprintf('%d generations completed\n',i);
%     end
%%
    
if M == 2 && Index ==2
    plot(chromosome(:,V + 1),chromosome(:,V + 2),'*');
    hold on
elseif M ==3 && Index ==2
    plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'*');
end    


if M == 2 && Index ==3
    subplot(1,2,1)
    plot([chromosome(:,V + 1)' range],[chromosome(:,V + 2)' y],'.');
    box off
    if i==1
    ax=max(chromosome(:,V + 1));
    ay=max(chromosome(:,V + 2));
    bx=min([chromosome(:,V + 1)' range]);
    by=min([chromosome(:,V + 2)' y]);
    end
    axis([bx ax by ay])
    title([num2str(i),' generations'])
    subplot(1,2,2)
    plot(1:size(IGD,2),IGD)
    title(['IGD'])
    box off
    pause(dt);
    drawnow;
elseif M ==3 && Index ==3
    plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'.');
    pause(dt);
    drawnow;
end    


if M == 2 && Index ==4
    subplot(1,2,1)
    plot([chromosome(:,V + 1)' range],[chromosome(:,V + 2)' y],'.');
    box off
    if i==1
    ax=max(chromosome(:,V + 1));
    ay=max(chromosome(:,V + 2));
    bx=min([chromosome(:,V + 1)' range]);
    by=min([chromosome(:,V + 2)' y]);
    end
    axis([bx ax by ay])
    title([num2str(i),' generations'])
    subplot(1,2,2)
    plot(1:size(IGD,2),IGD)
    title(['IGD'])
    box off
    pause(dt);
    drawnow;
    F0=getframe(gcf);
    I0=frame2im(F0);
    [I0,map]=rgb2ind(I0,256);
    imagename=['nsgaII_',name,'_',num2str(gen),'gen_',num2str(pop),'pop.gif'];
    if pic_num == 1
        imwrite(I0,map,imagename,'gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I0,map,imagename,'gif','WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;   
elseif M ==3 && Index ==4
    plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'.');
    pause(dt);
    drawnow;
end    
   
end

%% Result
% Save the result in ASCII text format.
save solution.txt chromosome -ASCII

if M == 2 && Index ==1
    plot(chromosome(:,V + 1),chromosome(:,V + 2),'*');
elseif M ==3 && Index ==1
    plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'*');
end    
%% Visualize
% The following is used to visualize the result if objective space
% dimension is visualizable.
if M == 2 && Index ==1
    subplot(1,2,1)
    plot(chromosome(:,V + 1),chromosome(:,V + 2),'*');
    title([num2str(gen),' generations'])
    subplot(1,2,2)
    plot(1:size(IGD,2),IGD);
    title(['IGD'])
elseif M ==3 && Index ==1
    plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'*');
end

fileID1 = fopen('judgement.txt','a+');
fprintf(fileID1,'%s%d%s%d%s%s%s%d%s%2.8f\n%s' , 'pop=',pop,' gen=',gen,' function:',name,'ii=',gen,'IGD=',IGD(gen),'                 ');
fprintf(fileID1,'%d                ',judgementrecord);
fprintf(fileID1,'\n%s','IGD:   ');
fprintf(fileID1,'%2.8f ',IGD);
fprintf(fileID1,'\n\n');
fclose(fileID1);    

end     
 

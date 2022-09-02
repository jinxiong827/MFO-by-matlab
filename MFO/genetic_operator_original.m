function f  = genetic_operator_original(parent_chromosome, M, V, mu, mum, l_limit, u_limit)
%mu=20,mum=20
[N1,m] = size(parent_chromosome);%N�ǽ�����еĸ�������
 
clear m
p = 1;
was_crossover = 0;%�Ƿ񽻲��־λ
was_mutation = 0;%�Ƿ�����־λ
 
for i = 1 : N1%������Ȼѭ��N�Σ�����ÿ��ѭ�������и��ʲ���2������1���Ӵ����������ղ������Ӵ�����������Լ��2N��
    if rand(1) < 0.9%�������0.9
        child_1 = [];
        child_2 = [];
        parent_1 =  max([unidrnd(N1) unidrnd(N1)]);
        parent_2 =  max([unidrnd(N1) unidrnd(N1)]);
        while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))
            parent_2 = unidrnd(N1);
        end
        parent_1 = parent_chromosome(parent_1,:);
        parent_2 = parent_chromosome(parent_2,:);
        for j = 1 : V
            u(j) = rand(1);
            if u(j) <= 0.5
                bq(j) = (2*u(j))^(1/(mu+1));
            else
                bq(j) = (1/(2*(1 - u(j))))^(1/(mu+1));
            end
            child_1(j) = 0.5*(((1 + bq(j))*parent_1(j)) + (1 - bq(j))*parent_2(j));
            child_2(j) = 0.5*(((1 - bq(j))*parent_1(j)) + (1 + bq(j))*parent_2(j));
            if child_1(j) > u_limit(j)
                child_1(j) = u_limit(j);
            elseif child_1(j) < l_limit(j)
                child_1(j) = l_limit(j);
            end
            if child_2(j) > u_limit(j)
                child_2(j) = u_limit(j);
            elseif child_2(j) < l_limit(j)
                child_2(j) = l_limit(j);
            end
        end
        was_crossover = 1;
        was_mutation = 0;
    else%if >0.9
        parent_3 =  unidrnd(N1);
        child_3 = parent_chromosome(parent_3,:);
        for j = 1 : V
           r(j) = rand(1);
           if r(j) < 0.5
               delta(j) = (2*r(j))^(1/(mum+1)) - 1;
           else
               delta(j) = 1 - (2*(1 - r(j)))^(1/(mum+1));
           end
           child_3(j) = child_3(j) + delta(j);
           if child_3(j) > u_limit(j) % ����Լ��
               child_3(j) = u_limit(j);
           elseif child_3(j) < l_limit(j)
               child_3(j) = l_limit(j);
           end
        end 
        was_mutation = 1;
        was_crossover = 0;
    end
    if was_crossover
        child(p,:) = child_1;
        child(p+1,:) = child_2;
        was_cossover = 0;
        p = p + 2;
    elseif was_mutation
        child(p,:) = child_3(1,1 : V );
        was_mutation = 0;
        p = p + 1;
    end
end
 
f = child;
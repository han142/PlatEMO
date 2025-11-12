function [ps,pf]=MOEADBDN(func_name,VRmin,VRmax,n_obj,Xpop,Max_Gen)


%%
%-----------------------------------------------------------------------------------------
%初始化
%-----------------------------------------------------------------------------------------
    %算法参数定义
    Nvar          = size(VRmin,2);      %决策变量数目
    Nobj          = n_obj;      % 目标数目
    Generations = Max_Gen;
    M = Nobj;
    if M == 2
        Xpop=100;
        H = Xpop-1;
    elseif M == 3
        Xpop=406;
        H = 27;
    elseif M == 4
        Xpop=560;
        H = 13;
    end
    %初始化权重向量
    [weight_size,weight_vector] = EqualWeight(H,M);%产生均匀权重  W:(Xpop*3)
    weight_vector(weight_vector==0) = 0.000001;
    %档案集合
    best_xreal = [];
    best_obj = [];
    b_size = 0;
    Boundary = [VRmax;VRmin];
    update_frequency = 0.05;                %权重更新频率
    freeze_time = 0.9;                      %当权重不允许改变时，进化的时间(比率)
    T = ceil(Xpop/10);                     %邻居数目
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
    % 初始随机种群
    Parent = [];  % 父代种群.
    Mutant = [];  % 突变种群.
    Child  = [];  % 子代种群.
    VRmin_1=repmat(VRmin,Xpop,1);
    VRmax_1=repmat(VRmax,Xpop,1);
    Parent=VRmin_1+(VRmax_1-VRmin_1).*rand(Xpop,Nvar); %初始化个体的位置
    for parent=1:Xpop
        JxParent(parent,1:Nobj) = feval(func_name,Parent(parent,:));
    end
    subParent1star=[Parent,JxParent];
    subParent1star(:,1:Nobj+Nvar+1) = non_domination_sort_mod(subParent1star,Nobj,Nvar);
    firstfront_all=length(find(subParent1star(:,Nobj+Nvar+1)==1));
    best_xreal(1:firstfront_all,1:Nvar)=subParent1star(1:firstfront_all,1:Nvar);
    best_obj(1:firstfront_all,1:Nobj)=subParent1star(1:firstfront_all,Nvar+1:Nvar+Nobj);
    b_size = firstfront_all;
    ideal_point = min(JxParent);
    Child = Parent;
    % 检测每个权值的邻居
    B = pdist2(weight_vector,weight_vector);
    [~,B] = sort(B,2);
    B = B(:,1:T);

%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
%%
    % 进化过程
    for Gene = 2 : Generations
        
        %----------------------------------------------------------------------------------------
        % 计算 k1_closest
        k1_closest = round(10 - 8 * Gene / Generations);
        % 计算距离矩阵
        distance1 = squareform(pdist(Parent));
        distance1(1:Xpop+1:end) = inf;
        % 计算第 k1_closest 个最近距离的中值
        random = sort(distance1, 2);
        random = random(:, k1_closest);
        sorted_random = sort(random);
        if mod(Xpop, 2) == 0
            radius1 = mean(sorted_random(Xpop/2-1:Xpop/2));
        else
            radius1 = sorted_random((Xpop+1)/2);
        end
        % 计算拥挤度
        crowd_degree1 = ones(Xpop, 1);
        for i = 1:Xpop
            neighbors = distance1(i, :) < radius1 & distance1(i, :) ~= inf;
            if any(neighbors)
                crowd_degree1(i) = prod(distance1(i, neighbors) / radius1);
            end
        end
        crowd_degree1 = 1 - crowd_degree1;
        % 选择个体
        [~, index_1] = sort(crowd_degree1, 'ascend');
        x_select = Parent(index_1(1:T), :);
        %----------------------------------------------------------------------------------------

        %-----------------------------------------------------------------------------------------
        for i=1:Xpop
            Parent_Child = Parent;
            a = rand;
            b = Gene/Generations;
            c = rand;
            if a < b
                if c>0.5
                    ScalingFactor1=0.8+0.4*rand;
                    P = B(i,randperm(size(B,2)));
                    Mutant(i,:)=Gen(Parent_Child(i,:),Parent_Child(P(1),:),Parent_Child(P(2),:),Boundary,ScalingFactor1);
                else
                    ScalingFactor1=0.8+0.4*rand;
                    P = randperm(Xpop);
                    Mutant(i,:)=Gen(Parent_Child(i,:),Parent_Child(P(1),:),Parent_Child(P(2),:),Boundary,ScalingFactor1);
                end
            else
                if c>0.5
                    rev1 = randperm(T);
                    ScalingFactor1=0.8+0.4*rand;
                    P = randperm(Xpop);
                    Mutant(i,:)=Gen(x_select(rev1(1),:),Parent_Child(P(1),:),Parent_Child(P(2),:),Boundary,ScalingFactor1);
                else
                    rev1 = randperm(T);
                    ScalingFactor1=0.8+0.4*rand;
                    Mutant(i,:)=Gen(Parent_Child(i,:),x_select(rev1(1),:),Parent_Child(i,:),Boundary,ScalingFactor1);
                end
            end
            for nvar=1:Nvar
                if Mutant(i,nvar)<VRmin(nvar)
                    Mutant(i,nvar) = (VRmin(nvar)+Parent_Child(i,nvar))/2;
                elseif Mutant(i,nvar)>VRmax(nvar)
                    Mutant(i,nvar) = (VRmax(nvar)+Parent_Child(i,nvar))/2;
                end
            end
            % Crossover operator
            Child(i,:) = Mutant(i,:);
            JxChild(i,1:Nobj) = feval(func_name,Child(i,:));
            ideal_point = min(ideal_point,JxChild(i,:));
            %-----------------------------------------------------------------------------------------

            %-----------------------------------------------------------------------------------------
            %更新MOEA/D中相邻子问题(权向量)的解
            if a < b
                if rand>0.1
                    P = B(i,randperm(size(B,2)));
                else
                    P = randperm(Xpop);
                end
            else
                if rand>0.1
                    P = B(index_1(rev1(1)),randperm(size(B,2)));
                else
                    P = randperm(Xpop);
                end
            end
            if (Gene<=Generations*freeze_time)
                nr = 1;
            else
                nr = 1;
                P = randperm(Xpop);
            end
            index_2 = randperm(length(P));
            d = 0;
            for j = 1:length(P)
                if d>=nr
                    break;
                end
                index = P(index_2(j));
                f1 = max(abs(JxParent(index,:)-ideal_point)./weight_vector(index,:));
                f2 = max(abs(JxChild(i,:)-ideal_point)./weight_vector(index,:));
                if f2<f1
                    %更新权向量对应的个体
                    %移除之前权重向量对应的个体中的权重向量
                    Parent(index,:) = Child(i,:);
                    JxParent(index,:) = JxChild(i,:);
                    %处理之前权重向量
                    d = d+1;
                elseif f2==f1
                    sum1 = sum(JxParent(index,:));
                    sum2 = sum(JxChild(i,:));
                    if sum2<sum1
                        %更新权向量对应的个体
                        %移除之前权重向量对应的个体中的权重向量
                        Parent(index,:) = Child(i,:);
                        JxParent(index,:) = JxChild(i,:);
                        d = d+1;
                    end
                end
            end
            %-----------------------------------------------------------------------------------------

            %-----------------------------------------------------------------------------------------
            % 归档
            % 假设 best_xreal 和 best_obj 已预分配
            % 检查支配关系
            if (~mod(Gene, Generations * freeze_time))
                select_index = [];
                for j = 1:b_size
                    f_ideal = abs(best_obj(j, :) - ideal_point);
                    sum3 = sum(f_ideal, 2);
                    weight_vector_random = f_ideal ./ sum3;
                    dis1 = pdist2(best_obj(j, :),best_obj);
                    [dis1,index1] = sort(dis1,2);
                    % 计算标量值
                    % 计算 PBI 标量值
                    theta = 200; % 惩罚因子
                    % 计算 d1（投影距离）
                    d1_parent = sum((best_obj(index1(2), :) - ideal_point) .* weight_vector_random) / norm(weight_vector_random);
                    d1_child = sum((best_obj(j, :) - ideal_point) .* weight_vector_random) / norm(weight_vector_random);
                    % 计算 d2（垂直距离）
                    proj_parent = d1_parent * weight_vector_random / norm(weight_vector_random);
                    proj_child = d1_child * weight_vector_random / norm(weight_vector_random);
                    d2_parent = norm(best_obj(index1(2), :) - ideal_point - proj_parent);
                    d2_child = norm(best_obj(j, :) - ideal_point - proj_child);
                    % 计算 PBI 值
                    f1 = d1_parent + theta * d2_parent;
                    f2 = d1_child + theta * d2_child;
    
                    if f2>f1
                        select_index = [select_index,j];
                    elseif f2==f1
                        sum1 = sum(JxParent(index1(2),:));
                        sum2 = sum(JxChild(i,:));
                        if sum2>sum1
                            select_index = [select_index,j];
                        end
                    end 
                end
                [delete_index,~] = sort(select_index,2,"descend");
                for j = 1:length(delete_index)
                    best_xreal(delete_index(j),:) = [];
                    best_obj(delete_index(j),:) = [];
                end
                b_size = b_size - length(delete_index);
            end

            dominates_child = any(all(best_obj(1:b_size, :) <= JxChild(i, :), 2));
            child_dominates = find(all(JxChild(i, :) <= best_obj(1:b_size, :), 2));
            if (Gene <= Generations * freeze_time)
                if ~dominates_child
                    % 添加新解
                    b_size = b_size + 1;
                    best_xreal(b_size, :) = Child(i, :);
                    best_obj(b_size, :) = JxChild(i, :);
                    if ~isempty(child_dominates)
                        % 删除被支配的解
                        best_xreal(child_dominates, :) = [];
                        best_obj(child_dominates, :) = [];
                        b_size = b_size - length(child_dominates);
                    end
                end
            else
                if ~dominates_child
                    f_ideal = abs(JxChild(i, :) - ideal_point);
                    sum3 = sum(f_ideal, 2);
                    weight_vector_random = f_ideal ./ sum3;
                    dis2 = pdist2(JxChild(i, :),best_obj);
                    [dis2,index2] = sort(dis2,2);
                    % 计算标量值
                    % 计算 PBI 标量值
                    theta = 200; % 惩罚因子
                    % 计算 d1（投影距离）
                    d1_parent = sum((best_obj(index2(1), :) - ideal_point) .* weight_vector_random) / norm(weight_vector_random);
                    d1_child = sum((JxChild(i, :) - ideal_point) .* weight_vector_random) / norm(weight_vector_random);
                    % 计算 d2（垂直距离）
                    proj_parent = d1_parent * weight_vector_random / norm(weight_vector_random);
                    proj_child = d1_child * weight_vector_random / norm(weight_vector_random);
                    d2_parent = norm(best_obj(index2(1), :) - ideal_point - proj_parent);
                    d2_child = norm(JxChild(i, :) - ideal_point - proj_child);
                    % 计算 PBI 值
                    f1 = d1_parent + theta * d2_parent;
                    f2 = d1_child + theta * d2_child;
    
                    if f2<f1
                        % 添加新解
                        b_size = b_size + 1;
                        best_xreal(b_size, :) = Child(i, :);
                        best_obj(b_size, :) = JxChild(i, :);
                        if ~isempty(child_dominates)
                            % 删除被支配的解
                            best_xreal(child_dominates, :) = [];
                            best_obj(child_dominates, :) = [];
                            b_size = b_size - length(child_dominates);
                        end
                    elseif f2==f1
                        sum1 = sum(JxParent(index2(1),:));
                        sum2 = sum(JxChild(i,:));
                        if sum2<sum1
                            %添加新解
                            b_size = b_size + 1;
                            best_xreal(b_size, :) = Child(i, :);
                            best_obj(b_size, :) = JxChild(i, :);
                            if ~isempty(child_dominates)
                                %删除被支配的解
                                best_xreal(child_dominates, :) = [];
                                best_obj(child_dominates, :) = [];
                                b_size = b_size - length(child_dominates);
                            end
                        end
                    end 
                end
            end
            %-----------------------------------------------------------------------------------------
        end

        k_closest = round(6-4*Gene/Generations);
        %-----------------------------------------------------------------------------------------
        %归档维护
        if b_size>Xpop
            best_obj_min = min(best_obj(:,1:Nobj));
            best_obj_max = max(best_obj(:,1:Nobj));
            for i = 1:M
                if best_obj_min(i) == best_obj_max(i)
                    for j = 1:b_size
                        best_obj_norm(j,i) = 0;
                    end
                else
                    for j = 1:b_size
                        best_obj_norm(j,i) = (best_obj(j,i)-best_obj_min(i))/(best_obj_max(i)-best_obj_min(i));
                    end
                end
            end
            distance_obj = [];
            for i = 1:b_size
                for j = i+1:b_size
                    distance_obj(i,j) = norm(best_obj_norm(i,:)-best_obj_norm(j,:));
                    distance_obj(j,i) = distance_obj(i,j);
                end
                distance_obj(i,i) = inf;
            end
            %计算归档维护半径
            if b_size == 1
                radius_obj = 0;
            else
                %求所有个体的第k个最近距离的中值
                [random,~] = sort(distance_obj,2);
                random = random(:,k_closest);
                [random,~] = sort(random,1);
                if mod(b_size,2) == 0
                    radius_obj = (random(b_size/2-1,1) + random(b_size/2,1))/2;
                else
                    radius_obj = random((b_size-1)/2,1);
                end
            end
            %个体拥挤度的初始化
            crowd_degree_obj = ones(b_size,1);
            clear neight_obj;
            neight_obj.index = [];
            neight_obj = repmat(neight_obj,1,b_size);
            %寻找邻居并计算拥挤度
            for i = 1:b_size
                for j = i+1:b_size
                    if distance_obj(i,j)<radius_obj
                        crowd_degree_obj(i,:) = crowd_degree_obj(i,:)*(distance_obj(i,j)/radius_obj);
                        crowd_degree_obj(j,:) = crowd_degree_obj(j,:)*(distance_obj(i,j)/radius_obj);
                        neight_obj(i).index = [neight_obj(i).index,j];
                        neight_obj(j).index = [neight_obj(j).index,i];
                    end
                end
            end
            %计算个体拥挤度
            crowd_degree_obj = 1-crowd_degree_obj;
            [max_obj_crowding,~] = max(crowd_degree_obj);

            best_xreal_min = min(best_xreal(:,1:Nvar));
            best_xreal_max = max(best_xreal(:,1:Nvar));
            for i = 1:Nvar
                if best_xreal_min(i) == best_xreal_max(i)
                    for j = 1:b_size
                        best_xreal_norm(j,i) = 0;
                    end
                else
                    for j = 1:b_size
                        best_xreal_norm(j,i) = (best_xreal(j,i)-best_xreal_min(i))/(best_xreal_max(i)-best_xreal_min(i));
                    end
                end
            end
            distance_xreal = [];
            for i = 1:b_size
                for j = i+1:b_size
                    distance_xreal(i,j) = norm(best_xreal_norm(i,:)-best_xreal_norm(j,:));
                    distance_xreal(j,i) = distance_xreal(i,j);
                end
                distance_xreal(i,i) = inf;
            end
            %计算归档维护半径
            if b_size == 1
                radius_xreal = 0;
            else
                %求所有个体的第k个最近距离的中值
                [random,~] = sort(distance_xreal,2);
                random = random(:,k_closest);
                [random,~] = sort(random,1);
                if mod(b_size,2) == 0
                    radius_xreal = (random(b_size/2-1,1) + random(b_size/2,1))/2;
                else
                    radius_xreal = random((b_size-1)/2,1);
                end
            end
            %个体拥挤度的初始化
            crowd_degree_xreal = ones(b_size,1);
            clear neight_xreal;
            neight_xreal.index = [];
            neight_xreal = repmat(neight_xreal,1,b_size);
            %寻找邻居并计算拥挤度
            for i = 1:b_size
                for j = i+1:b_size
                    if distance_xreal(i,j)<radius_xreal
                        crowd_degree_xreal(i,:) = crowd_degree_xreal(i,:)*(distance_xreal(i,j)/radius_xreal);
                        crowd_degree_xreal(j,:) = crowd_degree_xreal(j,:)*(distance_xreal(i,j)/radius_xreal);
                        neight_xreal(i).index = [neight_xreal(i).index,j];
                        neight_xreal(j).index = [neight_xreal(j).index,i];
                    end
                end
            end
            %计算个体拥挤度
            crowd_degree_xreal = 1-crowd_degree_xreal;
            [max_xreal_crowding,~] = max(crowd_degree_xreal);

            max_crowding = max_obj_crowding + max_xreal_crowding;
            flag = 0;
            delete_index_sum = [];
            if max_crowding == 0    %这意味着所有剩余的个体彼此之间不是相邻的
                %在这种情况下，随机删除一些，直到归档大小减少到容量
                num = b_size - Xpop;
                delete_index_sum = randperm(b_size);
                delete_index_sum = delete_index_sum(1,1:num);
                [delete_index_sum,~] = sort(delete_index_sum,2,"descend");
                for i = 1:num
                    best_obj(delete_index_sum(i),:) = [];
                    best_xreal(delete_index_sum(i),:) = [];
                end
                b_size = Xpop;
            else
                while b_size > Xpop
                    if max(crowd_degree_obj + crowd_degree_xreal)==0
                        flag = 1;
                        break;
                    else
                        [~,delete_index] = max(crowd_degree_obj + crowd_degree_xreal);
                        neight_obj_random = neight_obj(delete_index).index;
                        neight_xreal_random = neight_xreal(delete_index).index;
                        for i = 1:size(neight_obj_random,2)
                            neight_obj(neight_obj_random(i)).index(find(neight_obj(neight_obj_random(i)).index==delete_index)) = [];
                            crowd_degree_obj(neight_obj_random(i),:) = 1-((-(crowd_degree_obj(neight_obj_random(i),:)-1))*radius_obj/distance_obj(neight_obj_random(i),delete_index));
                        end
                        for i = 1:size(neight_xreal_random,2)
                            neight_xreal(neight_xreal_random(i)).index(find(neight_xreal(neight_xreal_random(i)).index==delete_index)) = [];
                            crowd_degree_xreal(neight_xreal_random(i),:) = 1-((-(crowd_degree_xreal(neight_xreal_random(i),:)-1))*radius_xreal/distance_xreal(neight_xreal_random(i),delete_index));
                        end
                        crowd_degree_obj(delete_index,:) = 0;
                        crowd_degree_xreal(delete_index,:) = 0;
                        delete_index_sum = [delete_index_sum,delete_index];
                        b_size = b_size - 1;
                    end
                end
                [delete_index_sum,~] = sort(delete_index_sum,2,"descend");
                for i = 1:size(delete_index_sum,2)
                    best_obj(delete_index_sum(i),:) = [];
                    best_xreal(delete_index_sum(i),:) = [];
                end
            end
            if flag == 1
                %在这种情况下，随机删除一些，直到归档大小减少到容量
                num = b_size - Xpop;
                delete_index_sum = randperm(b_size);
                delete_index_sum = delete_index_sum(1,1:num);
                [delete_index_sum,~] = sort(delete_index_sum,2,"descend");
                for i = 1:num
                    best_obj(delete_index_sum(i),:) = [];
                    best_xreal(delete_index_sum(i),:) = [];
                end
                b_size = Xpop;
            end
        end
        %-----------------------------------------------------------------------------------------

        
        %-----------------------------------------------------------------------------------------
        k2_closest = round(10 - 8 * Gene / Generations);
        if (~mod(Gene, Generations * update_frequency)) && (Gene <= Generations * freeze_time)
            %将best与parent结合起来选取最好的权重向量
            best_parent_xreal = [Parent;best_xreal];
            best_parent_obj = [JxParent;best_obj];
            subParent1star=[best_parent_xreal,best_parent_obj];
            subParent1star(:,1:Nobj+Nvar+1) = non_domination_sort_mod(subParent1star,Nobj,Nvar);
            firstfront_1=length(find(subParent1star(:,Nobj+Nvar+1)<subParent1star(Xpop,Nobj+Nvar+1)));
            firstfront_all=length(find(subParent1star(:,Nobj+Nvar+1)<=subParent1star(Xpop,Nobj+Nvar+1)));
            best_parent_xreal=subParent1star(1:firstfront_all,1:Nvar);
            best_parent_obj = subParent1star(1:firstfront_all,Nvar+1:Nvar+Nobj);
            sum_size = firstfront_all;

        
            % 计算距离矩阵
            distance2 = squareform(pdist(best_parent_obj));
            distance2(1:sum_size+1:end) = inf;
        
            % 计算第 k2_closest 个最近距离的中值
            random = sort(distance2, 2);
            random = random(:, k2_closest);
            sorted_random = sort(random);
            if mod(sum_size, 2) == 0
                radius2 = mean(sorted_random(sum_size/2-1:sum_size/2));
            else
                radius2 = sorted_random((sum_size+1)/2);
            end
        
            % 计算拥挤度和邻居
            crowd_degree2 = ones(sum_size, 1);
            neight = sparse(sum_size, sum_size);
            mask = distance2 < radius2 & distance2 ~= inf;
            neight(mask) = 1;
            neight = neight + neight';
            [row, col] = find(mask);
            for k = 1:length(row)
                i = row(k);
                j = col(k);
                crowd_degree2(i) = crowd_degree2(i) * (distance2(i,j) / radius2);
                crowd_degree2(j) = crowd_degree2(j) * (distance2(i,j) / radius2);
            end
            crowd_degree2 = 1 - crowd_degree2;
        
            % 删除个体
            if firstfront_1 < Xpop
                crowd_degree2(1:firstfront_1) = 0;
            end
            delete_index_sum = zeros(1, sum_size - Xpop);
            idx_count = 0;
        
            if max(crowd_degree2) == 0
                num = sum_size - Xpop;
                delete_index_sum = randperm(sum_size, num);
                idx_count = num;
            else
                while sum_size > Xpop
                    [~, delete_index] = max(crowd_degree2);
                    idx_count = idx_count + 1;
                    delete_index_sum(idx_count) = delete_index;
                    neight_idx = find(neight(delete_index, :));
                    neight(delete_index, :) = 0;
                    neight(:, delete_index) = 0;
                    for j = neight_idx
                        crowd_degree2(j) = 1 - ((1 - crowd_degree2(j)) * radius2 / distance2(j, delete_index));
                    end
                    crowd_degree2(delete_index) = 0;
                    sum_size = sum_size - 1;
                end
            end
        
            % 一次性删除
            if idx_count > 0
                delete_index_sum = sort(delete_index_sum(1:idx_count), 'descend');
                best_parent_obj(delete_index_sum, :) = [];
                best_parent_xreal(delete_index_sum, :) = [];
                sum_size = Xpop;
            end
        
            % 更新权重向量
            f_ideal = abs(best_parent_obj - ideal_point);
            sum3 = sum(f_ideal, 2);
            weight_vector = f_ideal ./ sum3;
        
            % 检测邻居
            B = pdist2(weight_vector, weight_vector);
            [~, B] = sort(B, 2);
            B = B(:, 1:T);
        
            Parent = best_parent_xreal;
            JxParent = best_parent_obj;
        end
        %-----------------------------------------------------------------------------------------

        %--------------------------------------------------------------------------------------------
        %cla;
        %Optimization_DrawGraph(Parent);%把函数值的点画到三维坐标系里。
        %Optimization_DrawGraph(best_obj(:,1:2));%把函数值的点画到三维坐标系里。
        %pause(0.01)
        %-----------------------------------------------------------------------------------------
    end
    %--------------------------------------------------------------------------------------------
    pf = best_obj;
    ps = best_xreal;
    %--------------------------------------------------------------------------------------------
    %figure(2)
    %Optimization_DrawGraph(Parent);%把函数值的点画到三维坐标系里。
    %figure(3)
    %Optimization_DrawGraph(JxParent);%把函数值的点画到三维坐标系里。
    %figure(4)
    %Optimization_DrawGraph(best_xreal);%把函数值的点画到三维坐标系里。
    %figure(5)
    %Optimization_DrawGraph(best_obj);%把函数值的点画到三维坐标系里。
    %--------------------------------------------------------------------------------------------
end


function [CA] = CSCD_ArchiveUpdate(Combination,N)
  [F,Front] = oursort(Combination.objs);
  [~,index] = CSCD(Combination.objs,Combination.decs,Front,F,10);
  CA = Combination(index(1:N));
end
function [CSCD,index] = CSCD(PopObj,PopDec,Front,F,KK)
    [N,M] = size(PopObj);
    D = size(PopDec,2);
% Sort the population according to the front number 对非支配层进行排序
    [Front,index_of_fronts] = sort(Front);
    PopObj = PopObj(index_of_fronts,:);
    PopDec = PopDec(index_of_fronts,:);
    finalpop = [];
    for front = 1 : (length(F) - 1)
        % 当前非支配层的个体
        Loc = find(Front == front);
        current_Pop = [PopDec(Loc,:) PopObj(Loc,:) Front(Loc)'];
        K = ceil(length(Loc)/KK);
        if length(Loc) < K
            keyboard
        end
        [ind,~] = kmeans(current_Pop(:,1:D),K);
        pop = cell(max(ind),1);
        for i1=1:length(ind)
            pop{ind(i1)} = [pop{ind(i1)};current_Pop(i1,:)];
        end
        for i2= 1 : size(pop,1)
            y = pop{i2};
            sorted_based_on_objective = [];
            for i = 1 : D
                [sorted_based_on_objective, index_of_objectives] = ...
                    sort(y(:,i));
                sorted_based_on_objective = [];
                for j = 1 : length(index_of_objectives)
                    sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
                end
                f_max = ...
                    sorted_based_on_objective(length(index_of_objectives), i);
                f_min = sorted_based_on_objective(1,  i);

                if length(index_of_objectives)==1
                    y(index_of_objectives(1),M + D + 1 + i) = 1;  %If there is only one point in current front
                elseif length(index_of_objectives)==2
                    y(index_of_objectives(1),M + D + 1 + i) = 1;
                    y(index_of_objectives(length(index_of_objectives)),M + D + 1 + i)=1;
                else
                    y(index_of_objectives(length(index_of_objectives)),M + D + 1 + i)...
                        = 2*(sorted_based_on_objective(length(index_of_objectives), i)-...
                        sorted_based_on_objective(length(index_of_objectives) -1, i))/(f_max - f_min);
                    y(index_of_objectives(1),M + D + 1 + i)=2*(sorted_based_on_objective(2, i)-...
                        sorted_based_on_objective(1, i))/(f_max - f_min);
                end
                for j = 2 : length(index_of_objectives) - 1
                    next_obj  = sorted_based_on_objective(j + 1, i);
                    previous_obj  = sorted_based_on_objective(j - 1,i);
                    if (f_max - f_min == 0)
                        y(index_of_objectives(j),M + D + 1 + i) = 1;
                    else
                        y(index_of_objectives(j),M + D + 1 + i) = ...
                            (next_obj - previous_obj)/(f_max - f_min);
                    end
                end
            end
            %% Calculate distance in decision space
            crowd_dist_var = [];
            crowd_dist_var(:,1) = zeros(size(pop{i2},1),1);
            for ii = 1 : D
                crowd_dist_var(:,1) = crowd_dist_var(:,1) + y(:,M + D + 1 + ii);
            end
            crowd_dist_var=crowd_dist_var./D;
            %avg_crowd_dist_var_k=mean(crowd_dist_var);
            y(:,M+D+2)=crowd_dist_var;
            y = y(:,1 : M + D +2 );
            finalpop=[finalpop;y];
        end
    end
        sorted_based_on_front = [];
%% calculate the objective space distance
[~,index_of_fronts] = sort(finalpop(:,M + D + 1));
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = finalpop(index_of_fronts(i),:);
   % avg_crowd_dist_var(i,:) = avg_crowd_dist_var(index_of_fronts(i),:);
end

  current_index = 0;

%% CSCD

    for front = 1 : (length(F) - 1)
  
        crowd_dist_obj = 0;
        y = [];
        previous_index = current_index + 1;
        for i = 1 : length(F(front).f)
            y(i,:) = sorted_based_on_front(current_index + i,:);%put the front_th rank into y
        end
        current_index = current_index + i;
   % Sort each individual based on the objective
        sorted_based_on_objective = [];
        for i = D+1 : M+D
            [sorted_based_on_objective, index_of_objectives] = ...
                sort(y(:,i));
            sorted_based_on_objective = [];
            for j = 1 : length(index_of_objectives)
                sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
            end
            f_max = ...
                sorted_based_on_objective(length(index_of_objectives), i);
            f_min = sorted_based_on_objective(1,  i);

            if length(index_of_objectives)==1
                y(index_of_objectives(1),M + D + 1 + i) = 1;  %If there is only one point in current front
            elseif i>D
                % deal with boundary points in objective space
                % In minimization problem, set the largest distance to the low boundary points and the smallest distance to the up boundary points
                y(index_of_objectives(1),M+D + 1 + i) = 1;
                y(index_of_objectives(length(index_of_objectives)),M+D + 1 + i)=0;
            end
             for j = 2 : length(index_of_objectives) - 1
                next_obj  = sorted_based_on_objective(j + 1, i);
                previous_obj  = sorted_based_on_objective(j - 1,i);
                if (f_max - f_min == 0)
                    y(index_of_objectives(j),M+D + 1 + i) = 1;
                else
                    y(index_of_objectives(j),M+D + 1 + i) = ...
                         (next_obj - previous_obj)/(f_max - f_min);
                end
             end
        end
    %% Calculate distance in objective space
        crowd_dist_obj = [];
        crowd_dist_obj(:,1) = zeros(length(F(front).f),1);
        for i = 1 : M
            crowd_dist_obj(:,1) = crowd_dist_obj(:,1) + y(:,M+D + 1+D + i);
        end
        crowd_dist_obj=crowd_dist_obj./M;
        avg_crowd_dist_obj=mean(crowd_dist_obj);
    %% Calculate special crowding distance
        special_crowd_dist=zeros(length(F(front).f),1);
        avg_crowd_dist_var = mean(y(:,M+D+2));
        for i = 1 : length(F(front).f)
            crowd_dist_var(i) = y(i,M+D+2);
            if crowd_dist_obj(i)>avg_crowd_dist_obj||crowd_dist_var(i)>avg_crowd_dist_var
                special_crowd_dist(i)=max(crowd_dist_obj(i),crowd_dist_var(i)); % Eq. (6) in the paper
            else
                special_crowd_dist(i)=min(crowd_dist_obj(i),crowd_dist_var(i)); % Eq. (7) in the paper
            end
        end
        y(:,M+D + 3) = special_crowd_dist;
        y(:,M+D+4)=crowd_dist_obj;
        [~,index_sorted_based_crowddist]=sort(special_crowd_dist,'descend');%sort the particles in the same front according to CSCD
        y=y(index_sorted_based_crowddist,:);
        y = y(:,1 : M+D+4 );
        z(previous_index:current_index,:) = y;
    end
    f= z();
    [CSCD,index] =sort( f(:,M + D+3));
end
function [F,Front] = oursort(PopObj)
    [N,M] = size(PopObj);
% Initialize the front number to 1.
    front = 1;
    
% There is nothing to this assignment, used only to manipulate easily in
% MATLAB.
    F(front).f = [];
    individual = [];
    Front = zeros(1,N);
%% Non-Dominated sort. 

    for i = 1 : N
        % Number of individuals that dominate this individual
        individual(i).n = 0; 
        % Individuals which this individual dominate
        individual(i).p = [];
        for j = 1 : N
            dom_less = 0;
            dom_equal = 0;
            dom_more = 0;
            for k = 1 : M
                if (PopObj(i,k) < PopObj(j,k))
                    dom_less = dom_less + 1;
                elseif (PopObj(i,k) == PopObj(j,k))  
                    dom_equal = dom_equal + 1;
                else
                    dom_more = dom_more + 1;
                end
            end
            if dom_less == 0 && dom_equal ~= M
                individual(i).n = individual(i).n + 1;
            elseif dom_more == 0 && dom_equal ~= M
                individual(i).p = [individual(i).p j];
            end
        end   
        if individual(i).n == 0
            Front(i) = 1;
            F(front).f = [F(front).f i];
        end
    end
% Find the subsequent fronts
    while ~isempty(F(front).f)
       Q = [];
       for i = 1 : length(F(front).f)
           if ~isempty(individual(F(front).f(i)).p)
                for j = 1 : length(individual(F(front).f(i)).p)
                    individual(individual(F(front).f(i)).p(j)).n = ...
                        individual(individual(F(front).f(i)).p(j)).n - 1;
                    if individual(individual(F(front).f(i)).p(j)).n == 0
                        Front(individual(F(front).f(i)).p(j)) = ...
                            front + 1;
                        Q = [Q individual(F(front).f(i)).p(j)];
                    end
               end
           end
       end
       front =  front + 1;
       F(front).f = Q;
    end
end
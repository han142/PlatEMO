function f = non_domination_sort_mod(x, M, V)

[N, m] = size(x);%此处的x指的是上面的particle
clear m

% Initialize the front number to 1.
front = 1;

% There is nothing to this assignment, used only to manipulate easily in
% MATLAB.
F(front).f = [];
individual = [];

for i = 1 : N
    % Number of individuals that dominate this individual该个体被支配的个数
    individual(i).n = 0; 
    % Individuals which this individual dominate该个体支配的个数
    individual(i).p = [];
    for j =1 : N
     
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M
            if (x(i,V + k) < x(j,V + k))
                dom_less = dom_less + 1;
            elseif (x(i,V + k) == x(j,V + k))  
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M
            individual(i).n = individual(i).n + 1;%被支配个体数
        elseif dom_more == 0 && dom_equal ~= M
            individual(i).p = [individual(i).p j];%所支配个体数
        end
    end   
    if individual(i).n == 0
        x(i,M + V + 1) = 1;%没有被支配，n+n_obj那一列置1.
        F(front).f = [F(front).f i];
    end
end
% Find the subsequent fronts寻找子Pareto前沿
while ~isempty(F(front).f)
   Q = [];
   for i = 1 : length(F(front).f)
       if ~isempty(individual(F(front).f(i)).p)
        	for j = 1 : length(individual(F(front).f(i)).p)
            	individual(individual(F(front).f(i)).p(j)).n = ...
                	individual(individual(F(front).f(i)).p(j)).n - 1;
        	   	if individual(individual(F(front).f(i)).p(j)).n == 0
               		x(individual(F(front).f(i)).p(j),M + V + 1) = ...
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end

[temp,index_of_fronts] = sort(x(:,M + V + 1));
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
end
current_index = 0;

f = sorted_based_on_front;
end

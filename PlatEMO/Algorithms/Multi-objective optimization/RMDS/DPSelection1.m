function [Choose,zmin_obj,zmin_dec,zmax_obj,zmax_dec] = DPSelection1(Population,zmin_obj,zmin_dec,zmax_obj,zmax_dec,K,t)
% The environmental selection of MDMMOEA
    
    PopObj = Population.objs;
    PopDec = Population.decs;
    [N,M]  = size(PopObj);
    [~,DM] = size(PopDec);
    A = 1 : N;         % Set of individuals wating for selection
    S = [];            % Set of selected individuals
    S_n = [];          % Set of neighbors

    %% Normalization
    % objective space  在临界层当中，求得每个目标上的最大值和最小值进行标准化
    zmin_obj = min(zmin_obj,min(PopObj,[],1));
    zmax_obj = max(zmax_obj,max(PopObj,[],1));
    PopObj = (PopObj - repmat(zmin_obj,N,1))./repmat(zmax_obj-zmin_obj,N,1);
    
    % decision space   在临界层当中，利用在决策空间上的最大值和最小值进行标准化
    zmin_dec = min(zmin_dec,min(PopDec,[],1));
    zmax_dec = max(zmax_dec,max(PopDec,[],1));
    PopDec = (PopDec - repmat(zmin_dec,N,1))./repmat(zmax_dec-zmin_dec,N,1);    
    
    %% Convergence and distribution indicators of each solution
    % covergence indicator
    c_obj = sqrt(sum(PopObj.^2,2));
    
    % decision space
    d_dec_e = pdist2(PopDec,PopDec);%最后一层当中所有个体之间的距离
    logical_matrix_dec = d_dec_e < 0.00000000001;%
    d_dec_num = sum(logical_matrix_dec,2);
    d_dec_num_new = d_dec_num;
    d_dec_num_new(d_dec_num==0) = 1;
    sum_e_dec = sum(d_dec_e .* logical_matrix_dec,2);
    fitness = (sum_e_dec) ./ (d_dec_num_new );
    
    %% Select the boundary solutions
     [num,~]=size(PopObj);
    for m = 1 : M
        if ~isempty(A)
            d_axis = sqrt(sum(PopObj(:,[1:m-1,m+1:end]).^2,2)) + c_obj;
            [A,S,S_n] = SelectingBoundarySolutions(d_axis,A,S,S_n,d_dec_e);
        end
    end
  
    
    %% Maximizing Diversity Selection
    while ~isempty(A)
        [A,S,S_n] = MaximizingDiversitySelection(fitness,A,S,PopObj,PopDec,d_dec_num_new,K,S_n,d_dec_e);
    end
    
    %% Make the number of selected solutions be K
    while length(S) < K
        if isempty(A)
            A = S_n;
            S_n = [];
        end
        [A,S,S_n] = MaximizingDiversitySelection(fitness,A,S,PopObj,PopDec,d_dec_num_new,K,S_n,d_dec_e);
    end
    
    Choose = S(1:K);
    

end

function [A,S,S_n] = MaximizingDiversitySelection(fitness,A,S,PopObj,PopDec,d_dec_num_new,K,S_n,d_dec_e)
        % calculate the repulsion distance
        d_obj = pdist2(PopObj(A,:),PopObj(S,:));
        d_dec = pdist2(PopDec(A,:),PopDec(S,:));
        d = d_obj + d_dec;
        temp = min(d,[],2);
        Fitness = temp;
    [~,x] = max(Fitness);
    x_p = A(x);
    S = [S,x_p];
    A(x) = [];
    % weaken the selection priority of neighbors
    xj = d_dec_e(A,x_p) < 0;
    if ~isempty(xj)
        S_n  = [S_n,A(xj)];
        A(xj) = [];
    end
end

function [A,S,S_n] = SelectingBoundarySolutions(d_axis,A,S,S_n,d_dec_e)
% Select the solution having the least d_axis value
    Fitness = d_axis(A);
    [~,x]   = min(Fitness);
    x_p = A(x);
    S = [S,x_p];
    A(x) = [];
    % weaken the selection priority of neighbors
    xj = d_dec_e(A,x_p) < 0;
    if ~isempty(xj)
        S_n  = [S_n,A(xj)];
        A(xj) = [];
    end
end

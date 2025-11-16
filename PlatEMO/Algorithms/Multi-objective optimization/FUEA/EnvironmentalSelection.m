
function [Population,CrowdDis] = EnvironmentalSelection(Population,N)
% EnvironmentalSelection - 局部收敛性排序并保留N个个体
% 输入：
%   Population : 候选种群
%   N          : 需要保留的个体数
% 输出：
%   Population : 筛选后的种群
%   CrowdDis   : 个体的拥挤距离

    n = length(Population);
    %% 计算局部收敛性
    dist = pdist2(Population.decs,Population.decs);
    V = 0.2*prod(max(Population.decs)-min(Population.decs)).^(1./size(Population.decs,2));
    DominationX = zeros(n); % 每对解的支配关系
    for i = 1:n
        for j = i+1:n
            if dist(i,j) > V
                continue;
            end
            L1 = Population(i).objs < Population(j).objs;
            L2 = Population(i).objs > Population(j).objs;
            if all(L1 | (~L2))
                DominationX(i,j) = 0;
                DominationX(j,i) = 1;
            elseif all(L2 | (~L1))
                DominationX(i,j) = 1;
                DominationX(j,i) = 0;
            end
        end
    end
    LocalC = zeros(1,n);
    for i = 1:n
        tmp = dist(i,:);
        index = tmp < V;
        LocalC(i) = (sum(DominationX(i,index))) ./ sum(index);
    end
    dist = sort(pdist2(Population.decs,Population.decs));
    CrowdDis = sum(dist(1:3,:));
    [~,index] = sortrows([LocalC' -CrowdDis']);
    Population = Population(index);
    if length(Population) > N
        Population = Population(1:N);
    end
    CrowdDis = Crowding(Population.decs);
end


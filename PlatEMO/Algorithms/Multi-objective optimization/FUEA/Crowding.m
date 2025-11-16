
function CrowdDis = Crowding(Pop)
% Crowding - 计算决策空间中解的调和平均距离
% 输入：Pop 为决策变量矩阵
% 输出：CrowdDis 为每个个体的拥挤距离

    [N,~] = size(Pop);
    K     = N - 1;
    Z     = min(Pop,[],1);
    Zmax  = max(Pop,[],1);
    pop   = (Pop - repmat(Z,N,1)) ./ repmat(Zmax - Z,N,1);
    distance = pdist2(pop,pop);
    [value,~] = sort(distance,2);
    CrowdDis  = K ./ sum(1./value(:,2:N),2);
end


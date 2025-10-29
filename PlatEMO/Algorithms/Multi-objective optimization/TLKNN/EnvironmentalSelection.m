function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N)
% TLKNN算法使用的环境选择算子

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% 非支配排序
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;

    %% 计算拥挤距离
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);

    %% 根据拥挤距离保留最后一层的解
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    remain   = min(N - sum(Next), numel(Last));
    if remain > 0
        Next(Last(Rank(1:remain))) = true;
    end

    %% 生成下一代种群
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

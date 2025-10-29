function [SubPopulation,c] = Partition(Population,K)
%Partition函数用于划分子种群，采用Kmeans截断方法
PopDec = Population.decs;
SubPopulation = struct('Ref', [], 'domains', []);
[ind,c]=kmeans(PopDec,K);%K是分割的个数
% c:簇质心位置，以数值矩阵形式返回。C 是 k×p 矩阵，其中第 j 行是簇 j 的质心。
% ind ：返回的是原矩阵的索引
for i=1:K
    SubPopulation(i).Ref = c(i,:);
    SubPopulation(i).domains = Population(ind==i);
end
end
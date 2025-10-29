function [SubPopulation,c] = Partition(Population,K)
PopDec = Population.decs;
SubPopulation = struct('Ref', [], 'domains', []);
[ind,c]=kmeans(PopDec,K);
for i=1:K
    SubPopulation(i).Ref = c(i,:);
    SubPopulation(i).domains = Population(ind==i);
end
end
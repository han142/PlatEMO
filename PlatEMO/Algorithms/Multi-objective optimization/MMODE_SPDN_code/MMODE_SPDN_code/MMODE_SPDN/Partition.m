function [SubPopulation,c] = Partition(Population,K)
%Partition�������ڻ�������Ⱥ������Kmeans�ضϷ���
PopDec = Population.decs;
SubPopulation = struct('Ref', [], 'domains', []);
[ind,c]=kmeans(PopDec,K);%K�Ƿָ�ĸ���
% c:������λ�ã�����ֵ������ʽ���ء�C �� k��p �������е� j ���Ǵ� j �����ġ�
% ind �����ص���ԭ���������
for i=1:K
    SubPopulation(i).Ref = c(i,:);
    SubPopulation(i).domains = Population(ind==i);
end
end
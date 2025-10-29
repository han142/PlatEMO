function Offspring = LocalSearch(Pop,EliteSolution,plotPops,D,M,ratio)
N = length(Pop);
Offspring = [];
Refs = [];
for i=1:N
    Ref = Pop(i).Ref;
    Refs = [Refs;Ref];
end
if ratio<0.5
    [ParseSolution,EIndex]=FindParseSolution(plotPops.decs); 
    Offspring = [Offspring MGA(EliteSolution,ParseSolution,[],{1,20,1,20,3})]; 
    [ParseSolution,EIndex]=FindParseSolution(EliteSolution.objs);
    CrowdedIndex = EIndex;
    Offspring = [Offspring MGA([EliteSolution(CrowdedIndex),Refs],[],EliteSolution,{1,20,1,20,2})]; 
end
[~,SelectIndex]=FindParseSolution(Refs); 
Parseindex = find(SelectIndex);
for i = 1:length(Parseindex)
    Parent = Pop(i).domains;
    Offspring = [Offspring MGA(Parent,Pop(i).Ref,EliteSolution,{1,20,1,20,1})]; 
end
a=1;
end
function [ParseSolution,SelectIndex]=FindParseSolution(PopDec)
N = size(PopDec,1);
for i=1:N
    if i==1
        distance(i)=dist(PopDec(i,:),PopDec(i+1,:)');
    elseif i==N
        distance(i)= dist(PopDec(i-1,:),PopDec(i,:)');
    else
        distance(i)= dist(PopDec(i-1,:),PopDec(i,:)')+dist(PopDec(i,:),PopDec(i+1,:)');
    end
end
SelectIndex = distance>(mean(distance)+std(distance));
ParseSolution = PopDec(SelectIndex,:);
end
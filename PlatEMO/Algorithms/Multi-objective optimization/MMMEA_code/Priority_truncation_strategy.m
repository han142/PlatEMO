function [Population,FrontNo,Priority] = Priority_truncation_strategy(Population,N)
    %% Non-dominated sorting
    Priority = ones(1,N);
    [~,ia,~] = unique(Population.objs,'rows');
    Population = Population(ia);
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;
    [Del1,Priority1]= Truncation(Population(FrontNo == MaxFNo).decs,sum(FrontNo == MaxFNo)-(N-sum(Next)));
    [Del2,Priority2]  = Truncation(Population(FrontNo == MaxFNo).objs,sum(FrontNo == MaxFNo)-(N-sum(Next)));
    Loc = ~all([Priority1' Priority2'],2);
    P1Max = max(Priority1);
    P2Max = max(Priority2);
    Priority1(~Del1) = Priority1(~Del1) + P1Max *P1Max;
    Priority2(~Del2) = Priority2(~Del2) + P2Max *P2Max;
    Priority = Priority1 .* Priority2;
    Del = [Del1',Del2'];
    DelLoc = find(any(Del,2));
    [~,index] = sort(Priority(DelLoc));
    Temp = find(FrontNo == MaxFNo);
    Next(setdiff(Temp,Temp(DelLoc(index(1:sum(FrontNo == MaxFNo)-(N-sum(Next))))))) = true;
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
end
function [Del,Priority] = Truncation(PopObj,K)
% Select part of the solutions by truncation
    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    Priority = zeros(1,size(PopObj,1));
    i=1;
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,index] = sortrows(Temp);
        Del(Remain(index(1))) = true;
        Priority(Remain(index(1))) = i;
        i = i+1;
    end
end
function Archive_Diversity = Centroid_truncation_strategy(Archive_Diversity,N)
N_Ref = length(Archive_Diversity);
if N_Ref>N
    for i=1:N_Ref
        Ref(i,:) = Archive_Diversity(i).Ref;
    end
    Del = Truncation(Ref,N_Ref-N);
    Archive_Diversity = Archive_Diversity(~Del);
end
end
function Del = Truncation(Ref,K)
% Select part of the solutions by truncation
    %% Truncation
    Distance = pdist2(Ref,Ref);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(Ref,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end
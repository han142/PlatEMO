function [Population,zmin_obj,zmin_dec,zmax_obj,zmax_dec] = EnvironmentalSelection1(Population,zmin_obj,zmin_dec,zmax_obj,zmax_dec,K,t)
% The environmental selection of MDMMOEA

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,K);
    Next = FrontNo < MaxFNo;
    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    [Choose,zmin_obj,zmin_dec,zmax_obj,zmax_dec] = DPSelection1(Population(Last),zmin_obj,zmin_dec,zmax_obj,zmax_dec,K-sum(Next),t);
    
    % Population for next generation
    Next(Last(Choose)) = true;
    Population = Population(Next);
    

end
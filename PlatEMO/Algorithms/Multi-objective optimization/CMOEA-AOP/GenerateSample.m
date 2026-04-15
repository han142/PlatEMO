function [state1, action1, nextstate1, ter, reward1] = GenerateSample(Problem, action,LastPopulation,Population)
%------------------------------- Reference --------------------------------
% Shao S, Tian Y, Yang S, et al. Deep Reinforcement Learning-Assisted 
% Automated Operator Portfolio for Constrained Multi-Objective Optimization[J]. 
% IEEE Transactions on Emerging Topics in Computational Intelligence, 2026.
%--------------------------------------------------------------------------    
    LastPopulationHV = HV(LastPopulation,max([LastPopulation.objs;Population.objs]));
    PopulationHV = HV(Population,max([LastPopulation.objs;Population.objs]));
    reward = (PopulationHV - LastPopulationHV)/LastPopulationHV;
    if isnan(reward)
        reward = 0;
    end

    LastObjsVar = var(LastPopulation.objs);
    LastObjsCon = sum(LastPopulation.objs);
    LastCV = sum(max(Population.cons,0),'all');
    LastRatio = Problem.FE/Problem.maxFE;
    LastState = [LastObjsVar, LastObjsCon, LastCV, LastRatio];

    CurrentObjsVar = var(Population.objs);
    CurrentObjsCon = sum(Population.objs);
    CurrentCV = sum(max(Population.cons,0),'all');
    CurrentRatio = Problem.FE/Problem.maxFE;
    CurrentState = [CurrentObjsVar, CurrentObjsCon, CurrentCV, CurrentRatio];

    state1 = LastState;
    action1 = action;
    nextstate1 = CurrentState;
    ter = 0;
    reward1 = reward;

end
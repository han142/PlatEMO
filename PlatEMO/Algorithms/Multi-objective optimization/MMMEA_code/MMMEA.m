function MMMEA(Global)
% <algorithm> <M>
% This is a simple demo of MMMEA
%--------------------------------------------------------------------------------------------------------
% If you find this code useful in your work, please cite the 
% following paper "H. Peng, S. Zhang, L. Li, B. Qu, X. Yue, Z. Wu, Multi-strategy multi-modal 
%		multi-objective evolutionary algorithm using macro and micro archive sets, 
%		Information Sciences 663 (2024) 120301".
%--------------------------------------------------------------------------------------------------------
% This function is implmented by SixiangZhang
%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
% More information can visit Hu Peng's homepage: https://whuph.github.io/index.html
%--------------------------------------------------------------------------------------------------------
    %% Generate the reference points and random population
    K = ceil(Global.N*0.5);
    Population   = Global.Initialization(Global.N*5);
    [FrontNo,~] = NDSort(Population.objs,Global.N*5);
    Archive_Convengence = Population;
    [SubPopulation,~] = Partition(Archive_Convengence,K);
    Archive_Diversity = SubPopulation;
    plotPops =[];
    for i=1:length(Archive_Diversity)
        plotPop = Archive_Diversity(i).domains;
        plotPops = [plotPops plotPop];
    end
    %% Optimization
    while Global.NotTermination(Archive_Convengence)
        Offspring = LocalSearch(Archive_Diversity,Archive_Convengence,plotPops,Global.D,Global.M,Global.gen/Global.maxgen);
        Combination = [Offspring,Archive_Convengence,plotPops];
        %% Update MIA
        [Archive_Convengence,FrontNo] = Priority_truncation_strategy(Combination,Global.N); 
        %% Update MAA
        [SubPopulation,~] = Partition(Archive_Convengence,K);
        Archive_Diversity = Centroid_truncation_strategy([SubPopulation,Archive_Diversity],Global.N);
        plotPops =[];
        for i=1:length(Archive_Diversity)
            plotPop = Archive_Diversity(i).domains;
            plotPops = [plotPops plotPop];
        end
    end
end
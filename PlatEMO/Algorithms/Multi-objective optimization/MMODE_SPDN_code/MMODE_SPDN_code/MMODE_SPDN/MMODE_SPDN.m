function  MMODE_SPDN( Global )
% <algorithm> <M>
% This is a simple demo of MMODE_SPDN
%--------------------------------------------------------------------------------------------------------
% If you find this code useful in your work, please cite the 
% following paper "H. Peng, W. Xia, Z. Luo, C. Deng, H. Wang, Z. Wu, A multimodal multi-objective 
% differential evolution with series-parallel combination and dynamic neighbor strategy, 
% Information Sciences,2024".
%--------------------------------------------------------------------------------------------------------
% This function is implmented by Wenwen Xia.
%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
% More information can visit Hu Peng's homepage: https://whuph.github.io/index.html
%--------------------------------------------------------------------------------------------------------
%% Parameter setting
div = Global.ParameterSet(10);
OS_neighbourhood = 25;
K=10;
eps = Global.ParameterSet(0.3);
%% Generate random population
Population = Global.Initialization();
Population1 = Global.Initialization();
temp1=Population;
temp2=Population;
EvoState = Global.gen / Global.maxgen;
[FrontNo,~] = NDSort(Population.objs,Global.N);
BA = Population(FrontNo(1:Global.N));
[DA,CrowdDis2] = ArchiveUpdate(temp2,Global.N,eps,0);
[CrowdDis,CrowdDec,~] = CrowdingDistance(Population.decs,Population.objs,FrontNo);
[SubPopulation,~] = Partition(BA,K);
ADA = SubPopulation;
%% Optimization
 while Global.NotTermination(Population1)
     Offspring1 = LocalSearch(ADA,Global.D,Global.M);
     Combination = [BA,Offspring1];
     MatingPoolAA = Selection(temp1,OS_neighbourhood,CrowdDis,CrowdDec);
     AA=Dynamic_neighbor_based_mutation(temp1(MatingPoolAA(:,1)),temp1(MatingPoolAA(:,2)),temp1(MatingPoolAA(:,3)),temp1(MatingPoolAA(:,4)),{0.1, 0.5, 1, 15});
     BA = CSCD_ArchiveUpdate(Combination,Global.N);
     MatingPoolCA = MatingSelection(temp2.decs,div);
     Offspring3 = GA(temp2(MatingPoolCA));
     CA = GrEA([Population,Offspring3],Global.N,div);
     MatingPool2 = TournamentSelection(2,round(Global.N),-CrowdDis2);
     Offspring2  = GA(DA(MatingPool2));
     [DA,CrowdDis2]    = Ranking([DA,Offspring2,CA],Global.N,eps,EvoState);
     Offspring  = Operator(Population);
     temp3=[AA,BA,CA,DA];
     Population1 = EnvironmentalSelection2([temp3,Offspring,Population1],Global.N*1.5,'Normal',EvoState);
     Population = EnvironmentalSelection1([temp3,Offspring,Population1],Global.N);
 end


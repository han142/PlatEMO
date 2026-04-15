classdef CMOEAAOP < ALGORITHM
% <2026> <multi> <real/integer/label/binary/permutation> <constrained>
% Automated Operator Portfolio Based CMOEA

%------------------------------- Reference --------------------------------
% Shao S, Tian Y, Yang S, et al. Deep Reinforcement Learning-Assisted 
% Automated Operator Portfolio for Constrained Multi-Objective Optimization[J]. 
% IEEE Transactions on Emerging Topics in Computational Intelligence, 2026.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% The original implementation of this work used a neural network written in Pytorch, 
% and the network design and code were developed with reference to the DDPG tutorial 
% available at https://hrl.boyuai.com/chapter/2/ddpg%E7%AE%97%E6%B3%95. To facilitate 
% open-source release and future comparisons by other researchers, we provide this 
% MATLAB implementation, which may differ from the original paper in some details.

    methods
        function main(Algorithm,Problem)
            Population{1} = Problem.Initialization();
            Population{2} = Problem.Initialization();
            Fitness{1}    = CalFitness(Population{1}.objs,Population{1}.cons);
            Fitness{2}    = CalFitness(Population{2}.objs);
            transfer_state=0;
            cnt=0;

            state_dim = Problem.M*2 + 2;      % 状态维度
            action_dim = 3;     % 动作维度
            action_bound = 5;   % 动作范围
            minnumber = 1e-5;
            
            [state1, ~, ~, ~, ~] = GenerateSample(Problem, rand(1,action_dim), Population{1}, Population{1});
            agent = DDPG(state_dim, action_dim, action_bound);

            %% Optimization
            while Algorithm.NotTerminated(Population{1})
                
                action = agent.selectAction(state1, true);
                action = (action+action_bound) / sum(action+action_bound+minnumber);
                cnt =cnt+1;                
                LastPopulation1 = Population{1};
                if transfer_state == 0
                    for i = 1: 2
                        valOffspring{i} = OperatorConstrainedAOP(Problem, Population, randperm(Problem.N,Problem.N), action, i);
                    end
                    
                    for i = 1:2
                        if i == 1
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{1},valOffspring{1},valOffspring{2}],Problem.N,i);
                        else
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{2},valOffspring{2},valOffspring{1}],Problem.N,i);
                        end
                    end
                    
                    if Problem.FE/Problem.maxFE >=0.2
                        transfer_state = 1;
                    end
                    
                else                    
                    for i = 1: 2
                        MatingPool = TournamentSelection(2,Problem.N,Fitness{i});                        
                        valOffspring{i} = OperatorConstrainedAOP(Problem, Population, MatingPool, action, i);
                    end
                    [~,~,Next] = EnvironmentalSelection( [Population{2},valOffspring{2}],Problem.N,1);
                    succ_rate(1,cnt) =  (sum(Next(1:Problem.N))/100) - (sum(Next(Problem.N+1:end))/50);
                    
                    [~,~,Next] = EnvironmentalSelection( [Population{1},valOffspring{1}],Problem.N,2);
                    succ_rate(2,cnt) =  (sum(Next(1:Problem.N))/100) - (sum(Next(Problem.N+1:end))/50);
                    
                    for i = 1:2
                        if   succ_rate(i,cnt) >0
                            rand_number = randperm(Problem.N);
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{i},valOffspring{i},Population{2/i}(rand_number(1:Problem.N/2))],Problem.N,i);
                        else
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{i},valOffspring{i},valOffspring{2/i}],Problem.N,i);
                        end
                    end
                end

                [state1, action1, nextstate1, ter, reward1] = GenerateSample(Problem, action, LastPopulation1, Population{1});
                agent.store(state1, action1, reward1, nextstate1, ter);
                agent.train();

            end
        end
    end
end
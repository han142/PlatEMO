classdef FUEA < ALGORITHM
% <2024> <multi> <real/integer> <multimodal>
% Opposition-based HREA: 加入对立学习与自适应差分进化的层次排序算法
% eps --- 0.3 --- 局部帕累托前的质量参数

%------------------------------- Reference --------------------------------
% 基于 W. Li 等人的 HREA 算法，增加对立学习初始化与自适应父代选择
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            eps = Algorithm.ParameterSet(0.3);
            p0  = 0.9;                        % 初始使用档案的概率
            %% 随机初始化并进行对立学习
             %% Generate the reference points    
            % R = UniformPoint(NR,Problem.M);
            Population = Problem.Initialization();
            OppDec     = repmat(Problem.upper,size(Population.decs,1),1)+ ...
                         repmat(Problem.lower,size(Population.decs,1),1)-Population.decs;
            OppPop     = Problem.Evaluation(OppDec);
            [Population,~] = EnvironmentalSelection([Population,OppPop],Problem.N);       
            [Archive,CrowdDis2]    = ArchiveUpdate(Population,Problem.N,eps,0);
            %% 进化迭代
            while Algorithm.NotTerminated(Archive)
                p = p0*(1-Problem.FE/Problem.maxFE);   % 随迭代次数减小
                if rand < p
                    %% 从档案中选择父代
                    MatingPool = TournamentSelection(2,round(Problem.N),-CrowdDis2);
                    P1         = Archive(MatingPool);
                    P2         = Archive(randi(length(Archive),1,length(MatingPool)));
                    P3         = Archive(randi(length(Archive),1,length(MatingPool)));
                    Offspring  = OperatorDE(Problem,P1,P2,P3); % 差分进化算子
                else
                    %% 从当前种群选择父代
                    Offspring = FuzzySearch2(Problem,Population,0.05);

                end
                %% 环境选择与档案更新
                [Population,~] = EnvironmentalSelection([Population,Offspring],Problem.N);
                [Archive,CrowdDis2]    = ArchiveUpdate([Archive,Offspring],Problem.N,eps,Problem.FE/Problem.maxFE);
            end
        end
    end
end



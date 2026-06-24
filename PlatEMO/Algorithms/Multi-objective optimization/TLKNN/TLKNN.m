classdef TLKNN < ALGORITHM
% <2023> <multi> <real/integer> <multimodal>
% 基于迁移学习的KNN引导多峰多目标优化算法
% K --- 5 --- 聚类中心的数量

%------------------------------- Reference --------------------------------
% 此算法演示了如何通过K近邻分类与环境选择相结合，维护多峰多目标
% 问题中的多个帕累托集合。
%--------------------------------------------------------------------------

    methods
        function main(Algorithm, Problem)
            K = Algorithm.ParameterSet(5);

            %% 生成初始种群并建立存档
            Population = Problem.Initialization();
            if Problem.N < K
                K = Problem.N;                % 确保聚类数不超过种群规模
            end
            Label = kmeans(Population.decs, K);                               % 对决策空间进行K均值聚类
            Mdl   = fitcknn(Population.decs, Label, 'NumNeighbors', 3);       % 训练KNN分类器
            Archive = Population;                                            % 初始化存档

            %% 进化迭代
            while Algorithm.NotTerminated(Archive)
                stage = Problem.FE / Problem.maxFE;                           % 当前进化进度 0~1

                % 根据阶段决定参与锦标赛的父代集合
                if stage <= 2/3
                    Pool = Population;                                       % 前期仅使用当前种群
                else
                    Pool = [Population, Archive];                            % 后期引入存档加强利用
                end
                [FrontNo,~] = NDSort(Pool.objs,Pool.cons,length(Pool));
                ObjCD = CrowdingDistance(Pool.objs,FrontNo);
                DecCD = CrowdingDistance(Pool.decs,FrontNo);
                if stage <= 1/3
                    Crowd = DecCD;                                           % 早期强调决策空间多样性
                elseif stage <= 2/3
                    Crowd = 0.5*DecCD + 0.5*ObjCD;                           % 中期兼顾两种多样性
                else
                    Crowd = ObjCD;                                           % 后期关注目标空间多样性
                end
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -Crowd);
                Offspring  = OperatorGA(Problem, Pool(MatingPool));

                % 预测子代的聚类标签
                LabelOff = predict(Mdl, Offspring.decs);

                % 合并父代与子代
                Combined = [Population, Offspring];
                Label    = [Label; LabelOff];

                % 按簇执行环境选择，维护多个局部PS
                NewPop = [];
                for k = 1:K
                    idx = find(Label==k);                    % 找到第k簇
                    if ~isempty(idx)
                        SubPop = Combined(idx);              % 取出该簇中的个体
                        % 仅在子种群规模允许的情况下执行环境选择
                        quota  = min(ceil(Problem.N/K), length(SubPop));
                        [SubPop,~,~] = EnvironmentalSelection(SubPop, quota);
                        NewPop = [NewPop, SubPop];
                    end
                end
                % 若新种群不足，则从全部个体中补足
                if length(NewPop) < Problem.N
                    remain = min(Problem.N - length(NewPop), length(Combined));
                    [Extra,~,~] = EnvironmentalSelection(Combined, remain);
                    NewPop = [NewPop, Extra];
                else
                    NewPop = NewPop(1:Problem.N);
                end
                Population = NewPop;                        % 更新种群

                % 更新存档，保持多样性
                Archive = ArchiveUpdate([Archive,Population], Problem.N, stage);

                % 重新训练分类器以捕捉新的聚类结构
                Label = kmeans(Population.decs, K);
                Mdl   = fitcknn(Population.decs, Label, 'NumNeighbors', 3);
            end
        end
    end
end

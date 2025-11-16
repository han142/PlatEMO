
function Offspring = FuzzySearch1(Problem,Population,ratio)
%RadarSearch - 全维雷达式模糊搜索，用于中后期增强多样性
% 输入：
%   Problem   : 问题对象
%   Population: 当前种群
%   ratio     : 参与搜索的比例 (0~1)
% 输出：
%   Offspring : 
%


    if nargin < 3
        ratio = 0.05;                     % 默认仅对 5% 个体进行搜索
    end
    N    = length(Population);
    selN = max(1,round(ratio*N));         % 被选中的个体数量
    Sel  = Population(randperm(N,selN));
    Dec  = Sel.decs;
    R    = Problem.upper - Problem.lower; % 决策空间范围
    sigma = 0.01*R;                       % 基本扰动幅度

    OffDec = [];
    for i = 1:selN
        center = Dec(i,:);
        for k = 1:8                      % 生成八个方向
            dir  = randn(1,Problem.D);   % 随机方向
            dir  = dir./norm(dir);
            rlen = 0.05 + 0.1*rand;      % 随机半径比例 5%~15%
            step = rlen.*(R.*dir);       % 按范围缩放的步长
            point = center + step + sigma.*randn(1,Problem.D);
            OffDec = [OffDec; point];
        end
    end

    % 边界处理
    OffDec = max(min(OffDec,repmat(Problem.upper,size(OffDec,1),1)),...
                     repmat(Problem.lower,size(OffDec,1),1));

    Offspring = Problem.Evaluation(OffDec);
end


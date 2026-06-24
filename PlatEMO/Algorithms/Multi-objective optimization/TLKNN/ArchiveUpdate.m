function Archive = ArchiveUpdate(Candidates,N,stage)
% TLKNN 存档更新：根据进化阶段平衡前沿与多样性
% Candidates : 参与筛选的个体集合
% N          : 存档的最大容量
% stage      : 当前进化进度 (0~1)

    % 非支配排序
    [FrontNo,MaxFNo] = NDSort(Candidates.objs,Candidates.cons,N);
    Next = FrontNo < MaxFNo;

    % 计算目标空间和决策空间的拥挤距离
    ObjCD = CrowdingDistance(Candidates.objs,FrontNo);
    DecCD = CrowdingDistance(Candidates.decs,FrontNo);

    % 根据进化阶段选择不同的多样性指标
    if stage <= 1/3
        Crowd = DecCD;                                % 早期强调决策空间多样性
    elseif stage <= 2/3
        Crowd = 0.5*DecCD + 0.5*ObjCD;                % 中期兼顾两种多样性
    else
        Crowd = ObjCD;                                % 后期关注目标空间多样性
    end

    % 保留最后一层中拥挤距离较大的个体
    Last = find(FrontNo==MaxFNo);
    [~,Rank] = sort(Crowd(Last),'descend');
    remain = min(N - sum(Next), numel(Last));
    if remain > 0
        Next(Last(Rank(1:remain))) = true;
    end

    % 更新存档
    Archive = Candidates(Next);
end

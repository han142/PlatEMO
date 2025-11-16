
function [Population,CrowdDis] = ArchiveUpdate(Population,N,eps,st)
% ArchiveUpdate - 更新外部档案并保持局部帕累托前
% 输入：
%   Population : 当前种群
%   N          : 档案容量
%   eps        : 判断局部帕累托前的阈值
%   st         : 进化阶段（0~1）
% 输出：
%   Population : 更新后的档案种群
%   CrowdDis   : 个体的拥挤距离

    n = length(Population);
    if eps~=1 && st<0.5
       eps = 2*(1-eps)/(2*st+1)+2*eps-1;
    end

    %% 选择全局帕累托前
    [FrontNo,MaxFNo] = NDSort(Population.objs,n);
    next        = FrontNo==1;
    first_pf    = Population(next);
    new_pop     = first_pf;
    remain_pop  = Population(~next);
    V = 0.2*prod(max(Population.decs)-min(Population.decs)).^(1./size(Population.decs,2));

    while ~isempty(remain_pop)
        %% 删除过近的解以保持多样性
        dist  = min(pdist2(new_pop.decs,remain_pop.decs));
        index = dist<V;
        remain_pop(index) = [];
        if isempty(remain_pop)
            break;
        end
        %% 选择剩余的解
        [FrontNo,MaxFNo] = NDSort(remain_pop.objs,length(remain_pop));
        pick_pop = remain_pop(FrontNo==1);
        [nF,~]   = NDSort([pick_pop.objs .* (1-eps); first_pf.objs],length(pick_pop)+length(first_pf));
        nF       = nF(1:length(pick_pop));
        maxnF    = max(nF);
        if maxnF>1
            new_pop    = [new_pop pick_pop(nF==1)];
            remain_pop = remain_pop(FrontNo~=1);
            break;
        else
            new_pop    = [new_pop pick_pop];
            remain_pop = remain_pop(FrontNo~=1);
        end
    end
    Population = new_pop;

    %% 平衡每层帕累托前的数量
    if length(Population) > N
        awd_index = [];
        [FrontNo,MaxFNo] = NDSort(Population.objs,length(Population));
        new_pop   = [];
        n_sub_pop = ceil(N/MaxFNo);
        sel_pop   = [];
        tmp_pop   = [];
        for i = 1:MaxFNo
            pop = Population(FrontNo==i);
            if length(pop) < n_sub_pop
                sel_pop   = [sel_pop pop];
                awd_index = [awd_index (n_sub_pop-length(pop)).*ones(1,length(pop))];
            else
                tmp_pop = [tmp_pop pop];
            end
        end
        while length(tmp_pop) > N - length(sel_pop)
            dist = pdist2(tmp_pop.decs,tmp_pop.decs);
            dist = sort(dist);
            dist = sum(dist(1:3,:),1);
            [~,ind] = min(dist);
            tmp_pop(ind) = [];
        end
        awd_index = [awd_index zeros(1,length(tmp_pop))] + 1;
        Population = [sel_pop tmp_pop];
        CrowdDis   = Crowding(Population.decs);
        CrowdDis   = CrowdDis .* awd_index';
    else
        CrowdDis = Crowding(Population.decs);
    end
end



function MatingPool = Selection(Population,OS_neighbourhood,CrowdDis,CrowdDec)
    PopObj = Population.objs;
    PopDec = Population.decs;
    N = size(PopObj,1);
    Obj_Distance = pdist2(PopObj,PopObj);
   Dec_Distance = pdist2(PopDec,PopDec);
 % [FrontNo,~] = NDSort(PopObj,N);
  % [Obj_Distance ,Dec_Distance,~] = CrowdingDistance(PopDec,PopObj,FrontNo); 
    [~,Obj_index] = sort(Obj_Distance,2);
    [~,Dec_index] = sort(Dec_Distance,2);
    p = rand(1,N);
    First_index = find(p < 0.33); % 第一种情况种群索引
    Second_index = find(p >= 0.33 & p < 0.66); % 第二种情况种群索引
    Third_index = setdiff(1:N,[First_index Second_index]); % 第三种情况种群索引
    MatingPool = zeros(N,4);
    %% 第一种情况——随机选择
    Obj_Neighbor_index = Obj_index(First_index,1:OS_neighbourhood); % 第一种情况种群邻居索引
    MatingPool(First_index,:) =  Obj_Neighbor_index(randi(OS_neighbourhood,length(First_index),4)); % 随机选择5个
    %% 第二种情况——根据决策空间的欧式距离来选择
    Dec_Neighbor_index = Dec_index(Second_index,1:OS_neighbourhood); % 第二种情况种群邻居索引
    Dec_Neighbor_CrowdDec = CrowdDec(Dec_Neighbor_index); % 确定种群邻居中对应的拥挤度
    Temp = randi(OS_neighbourhood,length(Second_index),4);
    Temp_Dec_Neighbor_index = Dec_Neighbor_index(Temp);
    Temp_Dec_Neighbor_CrowdDec = Dec_Neighbor_CrowdDec(Temp);
    [~,index] = sort(Temp_Dec_Neighbor_CrowdDec,2,'descend'); % 对选择的4个的拥挤度进行降序排序
    MatingPool(Second_index,:) = [Temp_Dec_Neighbor_index(index(:,1)),Temp_Dec_Neighbor_index(index(:,2:end))];
    %% 第三种情况——根据目标空间的欧式距离来选择
    Obj_Neighbor_index = Obj_index(Third_index,1:OS_neighbourhood); % 第三种情况种群邻居索引
    Obj_Neighbor_CrowdDis = CrowdDis(Obj_Neighbor_index); % 确定种群邻居中对应的拥挤度
    Temp = randi(OS_neighbourhood,length(Third_index),4);
    Temp_Obj_Neighbor_index = Obj_Neighbor_index(Temp);
    Temp_Obj_Neighbor_CrowdDis = Obj_Neighbor_CrowdDis(Temp);
    [~,index] = sort(Temp_Obj_Neighbor_CrowdDis,2,'descend'); % 对选择的5个的拥挤度进行降序排序
    MatingPool(Third_index,:) = [Temp_Obj_Neighbor_index(index(:,1)),Temp_Obj_Neighbor_index(index(:,2:end))];
end
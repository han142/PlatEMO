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
    First_index = find(p < 0.33); % ��һ�������Ⱥ����
    Second_index = find(p >= 0.33 & p < 0.66); % �ڶ��������Ⱥ����
    Third_index = setdiff(1:N,[First_index Second_index]); % �����������Ⱥ����
    MatingPool = zeros(N,4);
    %% ��һ������������ѡ��
    Obj_Neighbor_index = Obj_index(First_index,1:OS_neighbourhood); % ��һ�������Ⱥ�ھ�����
    MatingPool(First_index,:) =  Obj_Neighbor_index(randi(OS_neighbourhood,length(First_index),4)); % ���ѡ��5��
    %% �ڶ�������������ݾ��߿ռ��ŷʽ������ѡ��
    Dec_Neighbor_index = Dec_index(Second_index,1:OS_neighbourhood); % �ڶ��������Ⱥ�ھ�����
    Dec_Neighbor_CrowdDec = CrowdDec(Dec_Neighbor_index); % ȷ����Ⱥ�ھ��ж�Ӧ��ӵ����
    Temp = randi(OS_neighbourhood,length(Second_index),4);
    Temp_Dec_Neighbor_index = Dec_Neighbor_index(Temp);
    Temp_Dec_Neighbor_CrowdDec = Dec_Neighbor_CrowdDec(Temp);
    [~,index] = sort(Temp_Dec_Neighbor_CrowdDec,2,'descend'); % ��ѡ���4����ӵ���Ƚ��н�������
    MatingPool(Second_index,:) = [Temp_Dec_Neighbor_index(index(:,1)),Temp_Dec_Neighbor_index(index(:,2:end))];
    %% �����������������Ŀ��ռ��ŷʽ������ѡ��
    Obj_Neighbor_index = Obj_index(Third_index,1:OS_neighbourhood); % �����������Ⱥ�ھ�����
    Obj_Neighbor_CrowdDis = CrowdDis(Obj_Neighbor_index); % ȷ����Ⱥ�ھ��ж�Ӧ��ӵ����
    Temp = randi(OS_neighbourhood,length(Third_index),4);
    Temp_Obj_Neighbor_index = Obj_Neighbor_index(Temp);
    Temp_Obj_Neighbor_CrowdDis = Obj_Neighbor_CrowdDis(Temp);
    [~,index] = sort(Temp_Obj_Neighbor_CrowdDis,2,'descend'); % ��ѡ���5����ӵ���Ƚ��н�������
    MatingPool(Third_index,:) = [Temp_Obj_Neighbor_index(index(:,1)),Temp_Obj_Neighbor_index(index(:,2:end))];
end
function [outputArg1,outputArg2] = TSSKEA(inputArg1,inputArg2)%两阶段+自适应降维+0.01稀疏度条纹采样
%% MOEA-PSL identify edge-network biomarkers
%   此处显示详细说明
CANCER=cell(3,1);
CANCER{1,1}='BRCA';
CANCER{2,1}='LUSC';
CANCER{3,1}='LUAD';
path_in='E:\matlab\edge_net_data\';
path_out='E:\matlab\result\KGPDNB1\';
for cn=1:length(CANCER)
    cancer=CANCER{cn,1};
    where=strcat(path_in,cancer,'_edge_net\*.mat');
    numf = dir(where);
    experiment_num=30;
    Problem.maxFE=30000;
    Problem.N=300;
    Problem.lower=0;
    Problem.upper=1;
    for fnum = 1 :112
        if cn==3
            
            fnum=fnum+49;
        end
        data1=load([path_in,cancer,'_node_net\DE_',cancer,'_sample_result_',num2str(fnum),'.mat']);
%         data1=load([path_in,cancer,'_node_net\DE_LUNG_sample_result_',num2str(fnum),'.mat']);
%         test_adjacency=data1.subnetwork_adjacency;
        Prior_information=data1.label;
        data2=load([path_in,cancer,'_edge_net\DE_',cancer,'_edge_network_',num2str(fnum),'.mat']);
        test_adjacency=data2.edge_net;
        gg=data2.gg;
        Problem.D=size(test_adjacency,1); %维度
        test_net=zeros(Problem.D,Problem.D);
        test_net(test_adjacency~=0)=1;

        %% 处理先验信息（边网络）,找到先验节点所在的所有边
        non_zero_indices = find(Prior_information);
        Prior_information = false(size(gg, 1), 1);
        for i = 1:size(gg, 1)
            if any(ismember(gg(i, :), non_zero_indices))
                Prior_information(i) = true;
            end
        end

        %% paired edge :g
        Cons=CONS(test_net);
        Cons(all(Cons==0,2),:)=[];%删除全零行
%         cons=Cons;%约束条件
        Cnum=size(Cons,1);
        g=cell(Cnum,1);
        for i=1:Cnum
            g{i}=find(Cons(i,:)==1);
        end
        clear Cons
        %% 计算每个维度的目标函数值
        Dimension=eye(Problem.D,Problem.D);
        CV_num=Calcons(Problem.D,Cnum,g,Dimension);
        D_score=abs(CV_num-Cnum);
        Non_dominated_sol=cell(30,1);
        R=cell(30,1);
        for EXP_NUM=1:experiment_num
            CalNum=0;
            sLower = 0.99;
            sUpper = 0.999;
            %% Population initialization
            pi=Prior_information';
            pv = tiedrank(-D_score)'; 
            [~,ind]=find(pi~=0);
            pv(ind)=1; 


            [~,~,Mask] = vssps(Problem.N,Problem.D,sLower,sUpper);
            sumRows = sum(Mask, 2);  
            Delete = sumRows < 3;  
            Mask = Mask(~Delete, :);  
            [Mask,~] = unique(Mask,'rows');
            if size(Mask,1)<Problem.N
                newsoultion=Problem.N-size(Mask,1);
                ns=creatpop(newsoultion,Problem.D,D_score,g);
                Mask=[Mask;ns];
            end
            Dec = ones(Problem.N,Problem.D);
            Population = Dec.*Mask;
            [Population,calnum]  = SOLUTION(Population,test_adjacency,g,D_score);
            CalNum=CalNum+calnum;
            
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection_adr(Population,Dec,Mask,Problem.N,0,0);
            sv=zeros(1,Problem.D);
            Last_temp_num=0;

            %% Optimization
            rho = 0.5;
            RHO=0.5;
            while(CalNum<Problem.maxFE)
                delta= CalNum/Problem.maxFE;
                Popmarks   = extract_d(Population);
                MatingPool = TournamentSelection_hamming(length(Population),length(Population),Popmarks,FrontNo,CrowdDis);
                
                %-------------update sv-----------%
                First_Mask=Mask(FrontNo==1,:);
                [temp_num,~]=size(First_Mask);
                temp_vote=sum(First_Mask,1);
                sv(1,:)=(Last_temp_num/(Last_temp_num+temp_num))*sv(1,:)+(temp_num/(Last_temp_num+temp_num))*(temp_vote/temp_num);
                Last_temp_num=temp_num;

                % Operator
                if (delta/0.618) < 0.618 %% 第一阶段pv
                    Site=false(1,length(Population));
                    [rbm,allOne,other] = deal([]);
                    [OffMask,poss_s_num] = Operator_adrpv(MatingPool,Site,rbm,allOne,other,D_score,test_adjacency,pv,rho);
                else %% 第二阶段pv+sv+cv
                    if rho<0.5
                        index=FrontNo<ceil(max(FrontNo)/2);
                    else
                        index=FrontNo==1;
                    end

                    [rbm,allOne,other,Site] = trainRBM(Mask(index,:),Dec(index,:),rho,Problem.N); %训练RBM
                    if rand < 1/3
                        [OffMask,poss_s_num] = Operator_adrpv(MatingPool,Site,rbm,allOne,other,D_score,test_adjacency,pv,rho);
                    elseif rand > 2/3
                        [OffMask,poss_s_num] = Operator_adrsv(MatingPool,Site,rbm,allOne,other,D_score,test_adjacency,sv,rho);
                    else
                        [OffMask,poss_s_num] = Operator_adr(MatingPool,Site,rbm,allOne,other,D_score,test_adjacency,rho);
                    end
                end
                OffMask =logical(OffMask);
                OffDec  =ones(size(OffMask,1),size(OffMask,2));
                Offspring = OffDec.*OffMask;
                [Offspring,calnum]  = SOLUTION(Offspring,test_adjacency,g,D_score);
                CalNum=CalNum+calnum;
                    
                [Population,Dec,Mask,FrontNo,CrowdDis,sRatio] = EnvironmentalSelection_adr([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,length(Population),poss_s_num);
                rho = (rho+sRatio)/2;
                RHO=[RHO;rho];
            end
            PopObj      = extract(Population);
            Pop         = extract_d(Population);
            [FrontNo,~] = NDSort(PopObj,size(PopObj,1));
            outputpop=Pop(FrontNo==1,:);
            Non_dominated_sol{EXP_NUM}=outputpop;
            R{EXP_NUM}=RHO;
        end
        
        Non_dominated_sol1=cell2mat(Non_dominated_sol);
        
        [pop,~,~]=unique(Non_dominated_sol1,'rows');%将得到的非支配解进行去重
        
        [functionvalue,~] = Calfunctionvalue(pop,test_adjacency);%计算去重后的非支配解的目标函数
        
        [FrontNo,~] = NDSort(functionvalue,size(functionvalue,1));%对这些解进行pareto等级排序
        
        POP=pop(FrontNo==1,:);%输出最终的非支配解对应个体
        
        FV=functionvalue(FrontNo==1,:);
        
        [FVV,mod_position]=sortrows(FV);
        
        Final_Pop=POP(mod_position,:);% 输出最终排序过的个体
        
        %% 统计多模态出现的次数
        FVV(:,2)=-FVV(:,2);
        A=tabulate(FVV(:,1));
        B= A(:,2)~=0;
        c=A(B,:);

        %% 储存数据
            
        
        filename=strcat(path_out,cancer,'\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_TSSKEA_NDSol.mat');
        save(filename,'FVV');
        filename1=strcat(path_out,cancer,'\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_TSSKEA_Num.mat');
        save(filename1,'c');
        filename2=strcat(path_out,cancer,'\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_TSSKEA_POP.mat');
        save(filename2,'Final_Pop');
        filename3=strcat(path_out,cancer,'\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_TSSKEA_boxchart.mat');
        save(filename3,'Non_dominated_sol');
        filename4=strcat(path_out,cancer,'\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_TSSKEA_pho.mat');
        save(filename4,'R');

    end
end
end

function pop=creatpop(popnum,D,D_score,g)
%% 选择两个一定相连的基因
pop=zeros(popnum,D);
gg=cell2mat(g);
%% 剩下的随机选择
for i=1:popnum
    mustselectnum_position=randperm(size(gg,1),1);  %% 选择两个一定相连的基因置为1
    mustselectnum=gg(mustselectnum_position,:);
    pop(i,mustselectnum)=1;
    D_score(mustselectnum)=0;
    canditateD=find(D_score~=0);
    
    min_D_score=D_score;
    gennum=round(0.9*rand*size(canditateD,1));
    for j=1:gennum
        canditateD=find(min_D_score~=0);   %% 使之前选过的基因不在参与选择
        variables=randperm(size(canditateD,1),2);
        if D_score(canditateD(variables(1)))>D_score(canditateD(variables(2)))
            pop(i,canditateD(variables(1)))=1;
            
            min_D_score(canditateD(variables(1)))=0;
        else
            pop(i,canditateD(variables(2)))=1;
            
            min_D_score(canditateD(variables(2)))=0;
        end
    end
end
end
function A_adjacent=CONS(test_Net)
[z1,z2]=find(triu(test_Net)~=0);
z=[z1,z2];
NNN=length(test_Net);

N1=NNN;
[N2,~]=size(z);
%calculate the adjacency matrix of bipartite graph
A_adjacent=zeros(N2,N1);
for i=1:N2
    
    A_adjacent(i,z(i,1))=1;
    A_adjacent(i,z(i,2))=1;
    
end
end
function CVvalue = Calcons(popnum,Cnum,g,pop)
cv=zeros(Cnum,1);
CVvalue=zeros(popnum,1);
for i=1:popnum
    ind=pop(i,:);
    for j=1:Cnum
        cv(j)=max(1-sum(ind(g{j})),0);
    end
    CVvalue(i)=sum(cv);
end
end
function [Population,calnum]  = SOLUTION(Mask,test_adjacency,g,D_score)
[functionvalue,calnum1] = Calfunctionvalue(Mask,test_adjacency);
[Mask,tt]=FIX2(Mask,functionvalue,g,size(test_adjacency,2),D_score,test_adjacency);
[functionvalue,calnum2]=Calfunctionvalue_afterfix(Mask,test_adjacency,tt,functionvalue);
calnum=calnum1+calnum2;
Population=cell(1,size(Mask,1));
for i=1:size(Population,2)
    Population{1,i}.dec=Mask(i,:);
    Population{1,i}.obj=functionvalue(i,:);
    Population{1,i}.con=0;
    Population{1,i}.add=[];
end
end
function [Functionvalue, calnum]= Calfunctionvalue(pop,test_adjacency)
Functionvalue=zeros(size(pop,1),2);
functionvalue=zeros(size(pop,1),4);
for i=1:size(pop,1)
    Functionvalue(i,1)=sum(pop(i,:));
    %计算内部矩阵
    a=find(pop(i,:)==1);
    matrix=test_adjacency(a,:);
    inmatrix=matrix(:,a);
    genin=nonzeros(inmatrix);
    
    functionvalue(i,1)=abs(mean(genin));
    functionvalue(i,2)=std(genin);
    
    %计算外部矩阵
    gen=nonzeros(matrix);
    genout=setdiff(gen,genin);
    
    functionvalue(i,3)=abs(mean(genout));
    functionvalue(i,4)=functionvalue(i,1)*functionvalue(i,2)/functionvalue(i,3);
    
end
calnum=i;
Functionvalue(:,2)= -round(functionvalue(:,4)*100)/100;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,1)=size(test_adjacency,2);
Functionvalue(tt,2)=0;
end
function [Functionvalue,calnum] = Calfunctionvalue_afterfix(pop,test_adjacency,tt,functionvalue)

for i=1:size(tt,1)
    functionvalue(tt(i),:)=Calfunctionvalue(pop(tt(i),:),test_adjacency);
end
calnum=i;
Functionvalue=functionvalue;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,2)=0;
ttt= Functionvalue(:,2)==0;
Functionvalue(ttt,1)=size(test_adjacency,2);
end
function Parents = TournamentSelection_hamming(K,N,Population,FrontNo,SpCrowdDis)
index=zeros(K,1);
% Parents=zeros(K,size(Population,2));
hamming_dist=pdist2(Population,Population,'hamming');
[~,site2]=sort(hamming_dist,2);
K1=randperm(N,K);  %第一个父母编号
K2=site2(K1,2);    %第二个父母编号

for i=1:K
    if FrontNo(K1(i))<FrontNo(K2(i))
        index(i)=K1(i);
    elseif FrontNo(K1(i))==FrontNo(K2(i))
        if SpCrowdDis(K1(i))>=SpCrowdDis(K2(i))
            index(i)=K1(i);
        else
            index(i)=K2(i);
        end
    else
        index(i)=K2(i);
    end
    
end
Parents=Population(index,:);
end
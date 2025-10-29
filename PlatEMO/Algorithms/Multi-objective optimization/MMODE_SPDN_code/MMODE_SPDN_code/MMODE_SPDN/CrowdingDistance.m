function [CrowdDis,CrowdDec,SCD] = CrowdingDistance(PopDec,PopObj,FrontNo)
% Calculate the crowding distance of each solution front by front

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    CrowdDec = zeros(1,N);
    SCD = zeros(1,N);
    k = fix(0.02*N);
    [~,index] = sortrows(PopObj);
    Dis_Dec = pdist2(PopDec,PopDec);
    Dis_Dec(logical(eye(N))) = inf;
    [~,Dec_index] = sort(Dis_Dec,2);
    Fmax = max(PopObj,[],1);
    Fmin = min(PopObj,[],1);
    CrowdDec(1) = sum((k:-1:1).*Dis_Dec(1,Dec_index(1,1:k)));
    CrowdDec(end) = sum((k:-1:1).*Dis_Dec(end,Dec_index(end,1:k)));
    %% ¼ÆËãÓµ¼·¶È
    for i=2:N-1
        for j=1:M
            CrowdDis(index(1)) = rand;
            CrowdDis(index(end)) = rand;
            CrowdDis(index(i)) = CrowdDis(index(i)) + abs(PopObj(index(i-1),j)-PopObj(index(i+1),j))/(Fmax(j)-Fmin(j));
        end
        CrowdDec(i) = sum((k:-1:1).*Dis_Dec(i,Dec_index(i,1:k)));
    end
    avg_CrowdDis = mean(CrowdDis);
    avg_CrowdDec = mean(CrowdDec);
    %% ¼ÆËãSCD
    for i=1:N
        if CrowdDec(i) > avg_CrowdDec || CrowdDis(i) > avg_CrowdDis
            SCD(i) = max(CrowdDec(i),CrowdDis(i)/FrontNo(i));
            SCD(i)=var(SCD(i));
        else
            SCD(i) = min(CrowdDec(i),CrowdDis(i));
            SCD(i)=var(SCD(i));
        end
    end
end
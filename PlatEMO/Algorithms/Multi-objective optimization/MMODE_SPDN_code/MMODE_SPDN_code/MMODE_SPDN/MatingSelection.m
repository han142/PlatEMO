function MatingPool = MatingSelection(PopDec,div)

    N = size(PopDec,1);

    %% Calculate the grid location of each solution
    fmax = max(PopDec,[],1);
    fmin= min(PopDec,[],1);
    lb   = fmin-(fmax-fmin)/2/div;%
    ub   = fmax+(fmax-fmin)/2/div;%
    d    = (ub-lb)/div;
    lb   = repmat(lb,N,1);%repmat函数的真正作用是扩展所得到的值，将值扩展成矩阵
    d    = repmat(d,N,1);
    GLoc = floor((PopDec-lb)./d); 
    GLoc(isnan(GLoc)) = 0;
    
    %% Calculate the GD value of each solution
    GD = zeros(N)+inf;
    for i = 1 : N-1
        for j = i+1 : N
            GD(i,j) = sum(abs(GLoc(i,:)-GLoc(j,:)));%abs是绝对值
            GD(j,i) = GD(i,j);
        end
    end
    
    %% Calculate the GCD value of each solution
    GD  = max(size(PopDec,2)-GD,0);
    GCD = sum(GD,2);
    
    %% Binary tournament selection
    Parents1   = randi(N,1,N);
    Parents2   = randi(N,1,N);
    Dominate   = any(PopDec(Parents1,:)<PopDec(Parents2,:),2) - any(PopDec(Parents1,:)>PopDec(Parents2,:),2);
    GDominate  = any(GLoc(Parents1,:)<GLoc(Parents2,:),2) - any(GLoc(Parents1,:)>GLoc(Parents2,:),2);
    MatingPool = [Parents1(Dominate==1 | GDominate==1),...
                  Parents2(Dominate==-1 | GDominate==-1),...
                  Parents1(Dominate==0 & GDominate==0 & GCD(Parents1)<=GCD(Parents2)),...
                  Parents2(Dominate==0 & GDominate==0 & GCD(Parents1)>GCD(Parents2))];
end
function [rbm,allOne,other,Site] = trainRBM(Mask,Dec,UpdateRatio,N)

    allZero = all(~Mask,1);
    allOne  = all(Mask,1);
    other   = ~allZero & ~allOne;

    Site1 = false(1,floor(N/2));
    Site2 = false(1,floor(N/2));
    Site = [Site1,Site2];
    
    rbm = [];
    
    if rand < UpdateRatio %随机数小于后代存活率才使用RBM
        
        Site1(randperm(floor(N/2),randi(ceil(N/2*0.1)))) = true; % 从N中随机选择10%个体使用RBM生成子代
        
        Site2 = Site1;

        Site = [Site1,Site2];

        K = ceil(sum(mean(abs(Mask(:,other).*Dec(:,other))>1e-6,1)));
        
        K = min(K,ceil(size(Mask,2)*0.2));
        
        %rbm = RBM(sum(other),K,8,1,0,0.4,0.15);
        rbm = RBM(sum(other),K,10,1,0,0.5,0.1);
        
        rbm.train(Mask(:,other));

    end
end


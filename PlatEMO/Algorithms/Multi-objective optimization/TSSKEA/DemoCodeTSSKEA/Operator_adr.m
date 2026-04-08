function [OffMask,poss_s_num] = Operator_adr(ParentMask,Site,rbm,allOne,other,D_score,test_adjacency,UpdateRatio)
    
    Parent1Mask = ParentMask(1:floor(end/2),:);
    Parent2Mask = ParentMask(floor(end/2)+1:end,:);
    
    if any(Site)
        
        poss    = rbm.reduce(ParentMask(:,other));
        OffTemp = Cross1(sum(Site),size(poss,2),poss,delta);
        OffTemp = rbm.recover(OffTemp);
        
        OffMaskT1 = false(size(OffTemp,1),size(Parent1Mask,2));
        OffMaskT1(:,other)  = OffTemp;
        OffMaskT1(:,allOne) = true;
        
    else
        OffMaskT1 = [];
    end
        OffMaskT2 = false(sum(~Site),size(ParentMask,2));
        OffMaskT2(:,allOne) = true;
        OffMaskT2(:,other) = Cross1(sum(~Site),size(ParentMask(~Site,other),2),ParentMask(~Site,other));
        OffMask = [OffMaskT1;OffMaskT2];
        [OffMask,poss_s_num]= Muation(OffMask,D_score,test_adjacency,Site);

end

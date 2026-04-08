function [OffMask,poss_s_num] = Operator_adrsv(ParentMask,Site,rbm,allOne,other,D_score,test_adjacency,sv,UpdateRatio)
    
    Parent1Mask = ParentMask(1:floor(end/2),:);
    Parent2Mask = ParentMask(floor(end/2)+1:end,:);
    Site1 = Site(:,1:floor(end/2));
    Site2 = Site(:,floor(end/2)+1:end);

    if any(Site)        
        poss    = rbm.reduce(ParentMask(:,other));
        OffTemp = Cross1(sum(Site),size(poss,2),poss);
        OffTemp = rbm.recover(OffTemp);
        
        OffMaskT1 = false(size(OffTemp,1),size(Parent1Mask,2));
        OffMaskT1(:,other)  = OffTemp;
        OffMaskT1(:,allOne) = true;
        
    else
        OffMaskT1 = [];
    end
        OffMaskT2 = false(sum(~Site)/2,size(ParentMask,2));
        OffMaskT2(:,allOne) = true;
        OffMaskT2(:,other) = Cross_sv(Parent1Mask(~Site1,other),Parent2Mask(~Site2,other),sv(1,other));
        OffMask = [OffMaskT1;OffMaskT2];
        [OffMask,poss_s_num]= Muation(OffMask,D_score,test_adjacency,Site);

end
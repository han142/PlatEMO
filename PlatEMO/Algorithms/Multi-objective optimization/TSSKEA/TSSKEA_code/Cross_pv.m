function [OffMask] = Cross_pv(Parent1Mask,Parent2Mask,pv)
    OffMask = Parent1Mask;
    for i = 1 : size(OffMask,1)
        if rand < 0.5
            index = find(Parent1Mask(i,:)&~Parent2Mask(i,:));
            index = index(TS(-pv(index)));
            OffMask(i,index) = 0;
        else
            index = find(~Parent1Mask(i,:)&Parent2Mask(i,:));
            index = index(TS(pv(index)));
            OffMask(i,index) = Parent2Mask(i,index);
        end
    end
end

function index = TS(pv)
% Binary tournament selection

    if isempty(pv)
        index = [];
    else
        index = TournamentSelection(2,1,pv);
    end
end
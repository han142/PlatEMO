function [OffMask] = Cross_sv(Parent1Mask,Parent2Mask,sv)
    OffMask = Parent1Mask;
    rate0   = sv;           % The probability that  0 inverts to  1
    rate1   = 1-rate0;      % The probability that  1 inverts to  0

    for i = 1 : size(OffMask,1)
        diff = find(Parent1Mask(i,:)~=Parent2Mask(i,:));
        temp_rate1=rate1(diff);
        temp_rate0=rate0(diff);
        rate = zeros(1,length(diff));
        rate(logical(OffMask(i,diff)))  = temp_rate1(logical(OffMask(i,diff)));
        rate(logical(~OffMask(i,diff))) = temp_rate0(logical(~OffMask(i,diff)));
        exchange  = rand(1,length(diff)) < rate;
        OffMask(i,diff(exchange))=~OffMask(i,diff(exchange)); 
    end
end
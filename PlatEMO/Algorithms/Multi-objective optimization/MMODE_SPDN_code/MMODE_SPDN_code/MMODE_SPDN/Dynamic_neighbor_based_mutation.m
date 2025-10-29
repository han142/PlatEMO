function Offspring = Dynamic_neighbor_based_mutation(Parent1,Parent2,Parent3,Parent4,Parameter)
    %% Parameter setting
    if nargin > 4
        [CR,F,proM,disM] = deal(Parameter{:});
    else
        [CR,F,proM,disM] = deal(1,0.5,1,20);
    end
    if isa(Parent1(1),'INDIVIDUAL')
        calObj  = true;
        Parent1 = Parent1.decs;
        Parent2 = Parent2.decs;
        Parent3 = Parent3.decs;
        Parent4 = Parent4.decs;
    else
        calObj = false;
    end
    [N,D]  = size(Parent1);
    Global = GLOBAL.GetObj();

    %% Differental evolution
    Site = rand(N,D) < CR;
    Offspring       = Parent1;
    Offspring(Site) = Parent2(Site) + F*(Parent3(Site)-Parent4(Site));

    %% Polynomial mutation
    Lower = repmat(Global.lower,N,1);
    Upper = repmat(Global.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    if calObj
        Offspring = INDIVIDUAL(Offspring);
    end
end
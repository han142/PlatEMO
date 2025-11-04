function MatingPool = DensityMatingSelection1(Population,zmin_dec,zmax_dec)
% The mating selection of MDMMOEA

    %% get the values of decision space
    PopDec = Population.decs;
    [N,~]  = size(PopDec);
    
    %% update zmin and zmax
    zmin_dec  = min(zmin_dec,min(PopDec,[],1));
    zmax_dec  = max(zmax_dec,max(PopDec,[],1));
    
    %% normalize the value in decision space
    PopDec = (PopDec - repmat(zmin_dec,N,1))./repmat(zmax_dec-zmin_dec,N,1);
    
    %% Calculate the density of each solution
    d_dec  = pdist2(PopDec,PopDec);
    d_dec  = sort(d_dec,2);
    d = d_dec;
    fitness_density = 1./(sum(d(:,2:ceil(end/10)),2)+1);

    %% Binary tournament selection
    MatingPool = TournamentSelection(2,N,fitness_density);
end
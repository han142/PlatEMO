classdef RMDS < ALGORITHM
% <multi> <real> <multimodal>
% the main framework of RMDS 

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            % Ideal points and Nadir points of two space
            zmin_obj = min(Population.objs,[],1);
            zmax_obj = max(Population.objs,[],1);
            zmin_dec = min(Population.decs,[],1);
            zmax_dec = max(Population.decs,[],1);
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = DensityMatingSelection1(Population,zmin_dec,zmax_dec);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,zmin_obj,zmin_dec,zmax_obj,zmax_dec] = EnvironmentalSelection1([Population,Offspring],zmin_obj,zmin_dec,zmax_obj,zmax_dec,Problem.N,Problem.FE/Problem.maxFE);
            end
        end
    end
end
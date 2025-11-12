classdef SOLUTION
    properties
        obj    % 目标函数值
        xreal  % 决策变量
    end
    
    methods
        function obj = SOLUTION(objs, xreals)
            if nargin > 0
                obj.obj = objs;
                obj.xreal = xreals;
            end
        end
    end
end
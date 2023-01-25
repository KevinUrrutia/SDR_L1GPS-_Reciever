classdef test_check < handle

    properties
        x = 0;
        y = 0;
        z = 0;
    end

    methods 
        function obj = plot_test(obj)
            disp('check_point')
            obj.x = 0:pi/100:2*pi;
            obj.y = sin(obj.x);
%             plot(obj.x,obj.y)
        end
    end
end


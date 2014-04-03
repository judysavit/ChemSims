classdef Simulators
    properties (Abstract)
        input;
        deltaInput;
        netModel;
    end
    methods (Abstract)
        [T,X] = simulate(netModel, tspan, X0);
    end
end

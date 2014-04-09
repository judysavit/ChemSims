classdef Simulators
    properties (Abstract)
    end
    methods (Abstract)
        [T,X] = simulate(netModel, tspan, X0);
    end
end

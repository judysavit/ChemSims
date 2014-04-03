classdef ODESimulator < Simulators
    properties (Abstract)
        input;
        deltaInput;
        netModel;
    end
    methods (Abstract)
        function sim = ODESimulator(netModel, input, deltaInput)
          sim.netModel = netModel;
          sim.input = input;
          sim.deltaInput = deltaInput;
        end
        function X0 = initialize(netModel, input)
            opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);
            [T,X] = ode45(netModel.ODE_fcn, [1:100], X0, opts, netModels.params, this.input);
            %TODO: ACTUALLY CHECK FOR CONVERGENCE TO SS - CHANGE TSPAN OR
            %ODE OPTIONS
        end
        function [T,X] = simulate(netModel, tspan, X0)
            opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);
            [T,X] = ode45(netModel.ODE_fcn, tspan, X0, opts, netModels.params, this.input-this.deltaInput);
            netModel.ode_Fcn
        end
         function [T,X] = simAndSave(netModel, tspan, X0)
            path = netModel.folder;
        end
    end
end
classdef ODESimulator < Simulators
    methods(Static)
        function sim = ODESimulator()
        end
        function X_ss = initialize(netModel, X0)
            opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);
            [T,X] = ode45(netModel.ODE_fcn, [1:100], X0, opts);
            %TODO: ACTUALLY CHECK FOR CONVERGENCE TO SS - CHANGE TSPAN OR
            %ODE OPTIONS
            X_ss = X(end,:);
        end
        function [T,X] = simulate( netModel, tspan, delInput, X_ss)
            netModel.input = netModel.input+delInput;
            opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);
            [T,X] = ode45(netModel.ODE_fcn, tspan, X_ss, opts);
            netModel.ode_Fcn
        end
    end
end
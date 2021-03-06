classdef NetModels < handle
    properties
        speciesList = {'A', 'B', 'C'};
        
        % simulation properties
        delayBeforeStim = 0.5; %sec
        input = -1;
        deltaInput = -1;
        X0 = [-1 -1 -1];
        initialized = false;
        simTimeDefault = 6;
    end
    methods
        function model = NetModels(params, input, deltaInput)
            if exist('params','var'); model.params = ternif(~any(params<0), params, -1);end
            if exist('input','var'); model.input = ternif(~(input<0), input, -1);end
            if exist('deltaInput','var'); model.deltaInput = ternif(~(deltaInput<0), deltaInput, -1);end
            model.initialized = false;
        end
        function [T,X]=initializeInput(this, input, InitialX)
            %[T,X] = initializeInput(this, input, InitialX)
            %   input and InitialX are optional arguments
            if exist('input','var')
                this.input = ternif(~(input<0), input, -1);
            else
                this.input = ternif(this.input<1, 1, this.input);
            end
            if exist('InitialX','var'); IC = InitialX; else IC = floor(this.params(end-2:end)./3); end;
            if any(IC<1); error('Initial X (before initialization) is less than 1'); end;
            opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);%, 'OutputFcn', this.ode_outputFcn);
            T = [];
            X = [];
            Eps = [1 1 1];
            while(any(Eps>0.1))
                [Ttemp,Xtemp] = ode15s(@(t,y) ode_fcn(this, t, y), [1 50], IC, opts);
                Eps = abs(Xtemp(end,:)-Xtemp(end-1,:));
                T = [T; Ttemp];
                X = [X; Xtemp];
                if length(T)>500
                    break;
                end
            end
            
            this.X0 = X(end,:);
            assert(length(this.X0) == length(this.speciesList));
            this.initialized = true;
        end
        
        function [T,X] = simulate(this, mode, deltaInput, simTime)
            %[T,X] = simulate(this, mode, deltaInput)
            %   options for mode are CME, langevin, ODE
            %   deltaInput is an optional argument, only use it if you want
            %       it to be different than the originally initialized
            %       deltaInput, which is an attribute of this class
            if ~this.initialized; disp('WARNING: Running simulation without initializing to Steady State'); end;
            if ~exist('mode','var'); error('MUST tell the simulator what simulation algorithm to use'); end;
            if ~ismember(mode, {'CME', 'langevin','ODE','ODE_afterSS'}); error('Illegal simulation mode'); end;
            if exist('deltaInput','var'); this.deltaInput= deltaInput; end;
            sim = ternif(exist('simTime','var'), simTime, this.simTimeDefault);
           
            assert(~any(this.X0==[-1 -1 -1]));
            switch mode
                case 'CME'
                    [T2,X2] = directMethod(this.stoich_matrix, @(X) prop_fcn(this, X), [0 sim], round(this.X0));
                case 'langevin'
                    [T2,X2] = langevinCustom(this.stoich_matrix, @(X) prop_fcn(this, X), [0 sim], this.X0);
                case {'ODE','ODE_afterSS'}
                    opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);%, 'OutputFcn',this.ode_outputFcn);
                    [T2,X2] = ode15s(@(t,y) ode_fcn(this, t, y), [0 sim], this.X0, opts);
                otherwise
                    error('Incorrect simulation mode');
            end
            T = T2;
            X = X2;
%             T = [-this.delayBeforeStim; T2];
%             X = [this.X0; X2];
        end
    end
end

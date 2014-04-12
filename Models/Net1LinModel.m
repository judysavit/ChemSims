classdef Net1LinModel < handle
    % Species:
    %   AA - active A
    %   BA - active B
    %   CA - active C
    % Reaction network (deterministic):
    %   dA/dt = k_1*u(t)(AT-AA)-k_2*AA
    %   dB/dt = k_3*AA*(BT-BA)-k_4*BA
    %   dC/dt = k_5*AA*(CT-CA)-k_6*BA*CA-k_7.*CA
    % Reaction network (stochastic, 7 rxns):
    %   c1: u + AI -> u + AA
    %   c2: AA -> AI
    %   c3: AA + BI -> AA + BA
    %   c4: BA -> BI
    %   c5: AA + CI -> AA + CA
    %   c6: BA + CA <-> BACA -> BA + CI
    %   c7: CA -> CI
    % Parameters:
    %   7 stochastic rate constants (c values) above
    %   3 total A, B, C constraints (AT, BT, CT)
    %
    % METHODS:
    %   model = Net1LinModel(paramSet)     
    %   initialize(this, input, deltaInput)
    %   [T,X] = simulate(this, mode, deltaInput)
    
    properties
        % network properties
        speciesList = {'A', 'B', 'C'};
        M = 7;
        stoich_matrix = [ 1 0 0 ; % prod. A
            -1 0 0 ; % deg. A
            0 1 0 ; % prod. B
            0 -1 0; % deg. B
            0 0 1 ; % prod. C
            0 0 -1 ; %inh. C by B
            0 0 -1]; %deg. C
        params = zeros(1,10);
        folder = 'Data Lin Impl BC Data';
        
        % simulation properties
        delayBeforeStim = 2; %sec
        input = -1;
        deltaInput = -1;
        X0 = [-1 -1 -1];
        initialized = false;
    end
    methods
        function model = Net1LinModel(params, input, deltaInput)
            if exist('params','var'); model.params = ternif(~any(params<0), params, -1);end
            if exist('input','var'); model.input = ternif(~(input<0), input, -1);end
            if exist('deltaInput','var'); model.deltaInput = ternif(~(deltaInput<0), deltaInput, -1);end
            model.initialized = false;
        end
        function changeParams(this, paramSet)
            if length(paramSet)~=10; error('Parameter vector has incorrect size for this model'); end
            this.params = paramSet; 
            this.initialized = false;
        end
        function initialize(this, input, deltaInput)
            if ~(this.initialized)
                if exist('input','var')
                    this.input = ternif(~(input<0), input, -1);
                else
                    this.input = 1;
                end
                if exist('deltaInput','var')
                    this.deltaInput = ternif(~(deltaInput<0), deltaInput, -1);
                else
                    this.deltaInput = 0.2;
                end
                
                IC = floor(this.params(8:10)./3);
                opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);%, 'OutputFcn', this.ode_outputFcn);
                [T,X] = ode15s(@(t,y) ode_fcn(this, t, y), [1 500], IC, opts);
                
                this.X0 = X(end,:);
                assert(length(this.X0) == length(this.speciesList));
                this.initialized = true;
            end
        end
        function [T,X] = simulate(this, mode, deltaInput)
            %[T,X] = simulate(this, mode, deltaInput)
            %   options for mode are CME, langevin, ODE
            %   deltaInput is an optional argument, only use it if you want
            %       it to be different than the originally initialized
            %       deltaInput, which is an attribute of this class
            if ~this.initialized; disp('WARNING: Running simulation without initializing to Steady State'); end;
            if ~ismember(mode, {'CME', 'langevin','ODE','ODE_afterSS'}); error('Illegal simulation mode'); end;
            
            if exist('deltaInput','var'); this.deltaInput = ternif(~(deltaInput<0), deltaInput, error('Negative deltaInput value')); end

            assert(~any(this.X0==[-1 -1 -1]));
            switch mode
                case 'CME'
                    [T, X] = directMethod(this.stoich_matrix, @(X) prop_fcn(this, X), [0 100], this.X0);
                case 'langevin'
                    [T,X] = langevinCustom(this.stoich_matrix, @(X) prop_fcn(this, X), [0 30], this.X0);%   Returns:
                case {'ODE','ODE_afterSS'}
                    opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);%, 'OutputFcn',this.ode_outputFcn);
                    [T,X] = ode15s(@(t,y) ode_fcn(this, t, y), [0 30], this.X0, opts);
                otherwise
                    error('Incorrect simulation mode');
            end
        end
        function a = prop_fcn(this, X)
            if this.input<0; error('No input has been initialized'); end
            p = this.params;
            u=this.input;
            c1 = p(1);
            c2= p(2);
            c3 = p(3);
            c4 = p(4);
            c5 = p(5);
            c6 = p(6);
            c7 = p(7);
            AT = p(8);
            BT = p(9);
            CT = p(10);
            
            A = X(1);
            B = X(2);
            C = X(3);
            a(1) = c1.*u.*(AT-A);
            a(2) = c2.*A;
            a(3) = c3.*A.*(BT-B);
            a(4) = c4.*B;
            a(5) = c5.*A.*(CT-C);
            a(6) = c6.*B.*C;
            a(7) = c7.*C;
            a = a';
        end
        function dXdt = ode_fcn(this, t, X)
            if this.input<0 || this.deltaInput<0; error('No input has been initialized'); end
            c = this.params;
            if this.initialized && t > this.delayBeforeStim;
                u = this.input+this.deltaInput;
            else
                u=this.input;
            end
            A = X(1);
            B = X(2);
            C = X(3);
            AT = c(8);
            BT = c(9);
            CT = c(10);
            dXdt = [c(1)*u.*(AT-A) - c(2).*A;
                c(3).*A.*(BT-B) - c(4).*B;
                c(5).*A.*(CT-C) - c(6).*B.*C - c(7).*C];
        end
        function status = ode_outputFcn(t, X, p, u)
            disp('ODE OUTPUT FCN DOESNT WORK YET');
            status = 0;
            %             if this.input<0; error('No input has been initialized'); end
            %             c = this.params;
            %             u=this.input;
            %             A = X(1);
            %             B = X(2);
            %             C = X(3);
            %             dXdt = [c(1)*u.*(AT-A) - c(2).*A;
            %                 c(3).*A.*(BT-B) - c(4).*B;
            %                 c(5).*A.*(CT-C) - c(6).*B.*C - c(7).*C];
            %             epsilon = 10^-5;
            %             maxdXdt = matmax(dXdt);
            %             status = (maxdXdt<epsilon);
        end
    end
end

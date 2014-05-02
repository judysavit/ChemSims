classdef Net1LinModel < NetModels
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
    %   model = Net1LinModel(params, input, deltaInput) ALL OPT INPUTS
    %   changeParams(this, paramSet)
    %   [T,X]=initializeInput(this, input)
    %   [T,X] = simulate(this, mode, deltaInput)
    
    properties
        % network properties
        M = 7;
        stoich_matrix = [ 1 0 0 ; % prod. A
            -1 0 0 ; % deg. A
            0 1 0 ; % prod. B
            0 -1 0; % deg. B
            0 0 1 ; % prod. C
            0 0 -1 ; %inh. C by B
            0 0 -1]; %deg. C
        params = zeros(1,10);
    end
    methods
        function model = Net1LinModel(params, input, deltaInput)
            model@NetModels(params, input, deltaInput);
        end
        function changeParams(this, paramSet)
            if length(paramSet)~=10; error('Parameter vector has incorrect size for this model'); end
            this.params = paramSet;
            this.initialized = false;
        end
      
        function a = prop_fcn(this, X)
            % init is a boolean that determines if this is a initialization
            % or simulation run. true for initialization.
            if (this.input<0 || this.deltaInput<0); error('No input has been initialized, cannot simulate without X0'); end

            u = this.input+this.deltaInput;
            
            p = this.params;
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
            
            A = ternif(X(1)<=0, 0, X(1));
            B = ternif(X(2)<=0, 0, X(2));
            C = ternif(X(3)<=0, 0, X(3));
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
            if this.input<0 || (this.initialized && this.deltaInput<0); error('No input has been initialized'); end
            if this.initialized; u = this.input+this.deltaInput; else u = this.input; end
            
            c = this.params;
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

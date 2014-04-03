classdef Net1LinModel_BC < NetModels
% Species:
%   AA - active A
%   BA - active B
%   CA - active C
%   BC - BC complex

% Reaction network (deterministic):
%   dA/dt = k_1*u(t)(AT-AA)-k_2*AA
%   dB/dt = k_3*AA*(BT-BA)-k_4*BA
%   dBC/dt = k5*B*C - (k6+k7)*BC
%   dC/dt = k_8*AA*(CT-CA)-k_5*B*C + k_6*BC - k_9*C

% Reaction network (stochastic, 7 rxns):
%   c1: u + AI -> u + AA
%   c2: AA -> AI
%   c3: AA + BI -> AA + BA
%   c4: BA -> BI
%   c5: AA + CI -> AA + CA
%   c6: BA + CA -> BACA
%   c7: BACA -> BA + CA
%   c8: BACA -> BA + CI
%   c9: CA -> CI
    properties 
        speciesList = {'A', 'B', 'BC', 'C'};
        M = 7;
        stoich_matrix = [ 1 0 0 0; % prod. A
            -1 0 0 0; % deg. A
            0 1 0 0; % prod. B
            0 -1 0 0; % deg. B
            0 0 0 1; % prod. C
            0 -1 1 -1; %BC goes into complex
            0 1 -1 1; %BC goes out of complex
            0 1 -1 0; %BC goes to inactive C product
            0 0 0 -1]; %C is degraded
        params = zeros(1,10);
        input = 1;
    end
    methods 
        function model = Net1LinModel_BC(paramSet,input)
            model.params = paramSet;
            model.input = input;
        end
        function a = prop_fcn(this, X)
            p = this.params;
            u = this.input;
            c1 = p(1);
            c2= p(2);
            c3 = p(3);
            c4 = p(4);
            c5 = p(5);
            c6 = p(6);
            c7 = p(7);
            c8 = p(8);
            c9 = p(9);
            AT = p(10);
            BT = p(11);
            CT = p(12);
            
            A = X(1);
            B = X(2);
            BC = X(3);
            C = X(4);
            a(1) = c1.*u.*(AT-A);
            a(2) = c2.*A;
            a(3) = c3.*A.*(BT-B);
            a(4) = c4.*B;
            a(5) = c5.*A.*(CT-C);
            a(6) = c6.*B.*C;
            a(7) = c7.*BC;
            a(8) = c8.*BC;
            a(9) = c9.*C;
            a = a';
        end
        function dXdt = ODE_fcn(this, t,X)
            p = this.params;
            u = this.input;
            c1 = p(1);
            c2= p(2);
            c3 = p(3);
            c4 = p(4);
            c5 = p(5);
            c6 = p(6);
            c7 = p(7);
            c8 = p(8);
            c9 = p(9);
            AT = p(10);
            BT = p(11);
            CT = p(12);
            
            A = X(1);
            B = X(2);
            BC = X(3);
            C = X(4);
            
            dXdt = [c1*u.*(AT-A) - c2.*A;
                c3.*A.*(BT-B) - c4.*B - c5.*B.*C + c6.*BC;
                 c5.*BC - (c6+c7).*BC;
                 c8.*A.*(CT-C)-c5.*BC + c6.*BC- c9.*C];
        end
    end
end

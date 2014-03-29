classdef Net1LinModel < NetModels
    properties 
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
        input = 1;
    end
    methods 
        function model = Net1LinModel(paramSet,input)
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
            AT = p(8);
            BT = p(9);
            CT = p(10);
            A = X(1);
            B = X(2);
            C = X(3);
            dXdt = [c1*u.*(AT-A) - c2.*A;
                c3.*A.*(BT-B) - c4.*B;
                c5.*A.*(CT-C) - c6.*B.*C - c7.*C];
        end
    end
end

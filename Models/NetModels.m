classdef NetModels
    properties (Abstract)
        M;
        stoich_matrix;
        speciesList;
    end
    methods (Abstract)
        a = prop_fcn(X);
        dXdt = ODE_fcn(t,X);
    end
end

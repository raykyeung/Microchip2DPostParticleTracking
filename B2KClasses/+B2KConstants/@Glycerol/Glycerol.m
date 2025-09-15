classdef Glycerol
    properties(Constant)
        % mu = 1.41; %[Pa*s or kg/(m*s^3))]
        % rho = 1.26; %[kg/m^3]
        % nu = B2KConstants.Glycerol.mu ./ B2KConstants.Glycerol.rho; %[N*m*s/kg or m^2/s]

        mu = DimVar(1.41,'Pa-s')
        rho = DimVar(1.26,'kg/m^3')
        nu = B2KConstants.Glycerol.mu ./ B2KConstants.Glycerol.rho
    end
end
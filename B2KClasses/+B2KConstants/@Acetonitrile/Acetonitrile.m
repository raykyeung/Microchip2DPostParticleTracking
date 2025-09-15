classdef Acetonitrile
    properties(Constant)
        % mu = 3.43*10^-4; %[Pa*s or (N/m^2)*s or (kg*m*s/s^2)/(m^2) or kg/(m*s))]
        % rho = 786; %[kg/m^3]
        % nu = B2KConstants.Acetonitrile.mu ./ B2KConstants.Acetonitrile.rho; %[N*m*s/kg or m^2/s]

        mu = DimVar(UC(3.43*10^-4,1*10^-6),'Pa-s');
        rho = DimVar(UC(786,1),'kg/m^3');
        nu = B2KConstants.Acetonitrile.mu ./ B2KConstants.Acetonitrile.rho;
        epsilon_r = UC(37.50,0.05);
    end
end
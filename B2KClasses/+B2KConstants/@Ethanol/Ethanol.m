classdef Ethanol
    properties(Constant)
        % mu = 1.08*10^-3; %[Pa*s or kg/(m*s^3))]
        % rho = 789; %[kg/m^3]
        % nu = B2KConstants.Ethanol.mu ./ B2KConstants.Ethanol.rho; %[N*m*s/kg or m^2/s]

        mu = DimVar(UC(1.08*10^-3,1*10^-5),'Pa-s');
        rho = DimVar(UC(789,1),'kg/m^3');
        nu = B2KConstants.Ethanol.mu ./ B2KConstants.Ethanol.rho;
        epsilon_r = UC(24.55,0.05);
    end
end
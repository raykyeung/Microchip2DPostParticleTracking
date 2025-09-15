classdef Methanol
    properties(Constant)
        % mu = 5.40*10^-4; %[Pa*s or kg/(m*s^3))]
        % rho = 762; %[kg/m^3]
        % nu = B2KConstants.Methanol.mu ./ B2KConstants.Methanol.rho; %[N*m*s/kg or m^2/s]

        mu = DimVar(UC(5.40*10^-4,1*10^-6),'Pa-s'); %Dynamic viscosity
        rho = DimVar(UC(792,1),'kg/m^3'); %Density
        nu = B2KConstants.Methanol.mu ./ B2KConstants.Methanol.rho; %Kinematic viscosity
        epsilon_r = UC(32.70,0.05); %Dielectric constant or relative permittivity
    end
end
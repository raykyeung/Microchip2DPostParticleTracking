classdef pChip
    properties(Constant)
        % H = 1*10^-4; %[m]
        % W = 6*10^-4; %[m]
        % equivVolSpherD = ((6 .* B2KConstants.pChip.W .* B2KConstants.pChip.W .* B2KConstants.pChip.H ./ pi) .^ (1/3))/1000; %[m]
        % D_max = (B2KConstants.pChip.W*sqrt(2))/1000; %[m]
        % sphericity = (((pi^(1/3))*(6*B2KConstants.pChip.W*B2KConstants.pChip.W*B2KConstants.pChip.H)^(2/3))/(4*B2KConstants.pChip.W*B2KConstants.pChip.H+2*B2KConstants.pChip.W*B2KConstants.pChip.W)) %[m]
        % density = 2330; %[kg/m^3] polysilicon
        % density = 2650; %[kg/m^3] silicon dioxide
        % density = 2329 %[kg/m^3] silicon

        H = DimVar(UC(1*10^-4,1*10^-6),'m');
        W = DimVar(UC(6*10^-4,1*10^-6),'m');
        EquivVolSpherD = (6 .* B2KConstants.pChip.W .* B2KConstants.pChip.W .* B2KConstants.pChip.H ./ pi) .^ (1/3);
        D_max = (B2KConstants.pChip.W .* sqrt(2))
        Sphericity = ((pi .^ (1/3)) .* (6 .* B2KConstants.pChip.W .* B2KConstants.pChip.W .* B2KConstants.pChip.H) .^ (2/3))...
            ./ (4 .* B2KConstants.pChip.W .* B2KConstants.pChip.H + 2 .* B2KConstants.pChip.W .* B2KConstants.pChip.W);
        Density = DimVar(UC(2329,1),'kg/m^3');
    end
end
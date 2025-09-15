classdef EdgePixel    
    properties
        position;       % 1D index inside image 
        x, y;           % subpixel position
        nx, ny;         % normal vector (normalized)
        curv;           % curvature
        i0, i1;         % intensities at both sides
    end

    % to know subscript of an edge pixel use command:
    % [row colum] = ind2sub(size(image),position)
    
    methods
    end 
end


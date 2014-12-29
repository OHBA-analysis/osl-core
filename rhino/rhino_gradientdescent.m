function M = rhino_gradientdescent(data1,data2)
% Registers two point clouds using gradient descent
% M = rhino_gradientdescent(data1,data2)
%
% REQUIRED INPUTS:
%
% data1      - [3 x N] locations of the set of static points 
%              (e.g. MRI scalp surface)
% data2      - [3 x N] locations of the set of moving points 
%              (e.g. Polhemus headshape points)
%
% OUTPUTS:
%
% M          - [4 x 4] rigid transformation matrix mapping data2 to data1
%
% Adam Baker 2014

    data1 = data1';
    data2 = data2';
    
    [x,fval] = fminsearch(@cost_function,zeros(1,6));
        
    M = spm_eeg_inv_rigidreg(data2r',data1');

    
    function err = cost_function(x)
        
        rx = x(1);
        ry = x(2);
        rz = x(3);
        tx = x(4);
        ty = x(5);
        tz = x(6);
    
        Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
        Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)];
        Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];
        T  = [tx ty tz];
        Mr = [Rx*Ry*Rz T'; 0 0 0 1];
                            
        data2r  = Mr*[data2; ones(1,size(data2,2))];
        data2r  = data2r(1:3,:);
                            
        err = sum(rhino_cperror(data1,data2r));
 
    end
    
end


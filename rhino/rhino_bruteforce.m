function Err = rhino_bruteforce(data1,data2)

% Adam Baker 2014
    
    maxRot   = 9; % degrees
    maxTrans = 10; % mm
  
    stepRot   = 3; % degrees
    stepTrans = 5; % mm
    
    gridRot   = (-maxRot:stepRot:maxRot)*pi/180;
    gridTrans = -maxTrans:stepTrans:maxTrans;

    data1 = data1';
    data2 = data2';
    
    m = size(data2,2);
    n = size(data1,2);
    
    Err = zeros(length(gridRot),length(gridRot),length(gridRot),length(gridTrans),length(gridTrans),length(gridTrans));
    
    for rx = gridRot
        for ry = gridRot
            for rz = gridRot
                for tx = gridTrans
                    for ty = gridTrans
                        for tz = gridTrans
                            
                            Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
                            Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)];
                            Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];
                            T  = [tx ty tz];
                            Mr = [Rx*Ry*Rz T'; 0 0 0 1];
                            
                            data2r  = Mr*[data2; ones(1,size(data2,2))];
                            data2r  = data2r(1:3,:);
                            
                            Err(gridRot==rx,gridRot==ry,gridRot==rz,gridTrans==tx,gridTrans==ty,gridTrans==tz) = sum(rhino_cperror(data1,data2r,'point'));
                            
                        end
                    end
                end
            end
        end
    end
    
end
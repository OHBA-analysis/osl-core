function [ret, w]=runcmd(cmd)

[ret, w]=system(cmd);

if(ret ~= 0),
    
    msg=['runcmd call: \n' cmd ' \nProduced error: \n' w];
    ME = MException('runcmd:error',msg);
    throw(ME);
          
end;

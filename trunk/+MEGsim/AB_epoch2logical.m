function l = AB_epoch2logical(ev,t)
% AB 2013

l = false(size(t));

for i=1:size(ev,1)
  l = l | ( t >= ev(i,1) & t < ev(i,2) );
end

end



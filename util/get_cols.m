function ret = get_cols( ind )

cols={'r','g','b','m','c','y','r--','g--','b--','m--','c--','y--'};

if nargin>0
    ret=cols{ind};
else
    ret=cols
end

end


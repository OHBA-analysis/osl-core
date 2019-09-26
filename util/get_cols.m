function ret = get_cols( ind )

    cols={'y--','r','g','b','m','c','y','r--','g--','b--','m--','c--'};

    if nargin>0

        ret=cols{rem(ind,length(cols))+1};

    else
        ret=cols;
    end

end


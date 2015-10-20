function pretty_string(str)
% Print a centered string with dashes above and below
strlen = 70;
padlen = (strlen - length(str))/2;
str = [repmat(sprintf(' '),1,floor(padlen)) str repmat(sprintf(' '),1,ceil(padlen))];
padstr = repmat(sprintf('-'),1,strlen);
fprintf('\n%s\n%s\n%s\n',padstr,str,padstr);

end


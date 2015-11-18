function x = read_vest(filename)

fp=fopen(filename,'r');

while(1)
	line = fgetl(fp);
	if(strcmp(line,'/Matrix'))
		break;
	end;
end;

str = fgetl(fp);
x = read_vest_line(str);
str = fgetl(fp);

c=1;
while(isstr(str))
	c=c+1;
	x(c,:) = read_vest_line(str);
	str = fgetl(fp);
end;

fclose(fp);

%%%%%%%%%%%%%%%%

function val = read_vest_line(str)

indDelim = sort([0,findstr(str,sprintf('\t')),findstr(str,sprintf(' '))]);
val = zeros(1,length(indDelim)-1);
s=0;

for ctDelim = 1:length(indDelim)-1,
      s=s+1;
      val(s) = str2double(str(indDelim(ctDelim)+1:indDelim(ctDelim+1)-1));
end

if any(isnan(val)),
   val=NaN;
end


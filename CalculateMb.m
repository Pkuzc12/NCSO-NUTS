function [Mb] = CalculateMb(x,size)

Mb=0;

for i=1:1:size
    for j=1:1:size
        index=((i-1)*size+j-1)*4;
        Mb=Mb+sin(x(index+1))*sin(x(index+2))+sin(x(index+3))*sin(x(index+4));
    end
end

Mb=Mb/(size*size*2);


end


function [A] = ChartGrt(size)


for x=1:1:size
    for y=1:1:size
        A(x,y).left.neighbor=[[IndexNext(x,size,-1),y];[x,IndexNext(y,size,-1)];[x,y]];
        A(x,y).right.neighbor=[[IndexNext(x,size,1),y];[x,IndexNext(y,size,1)];[x,y]];
    end
end


for x=1:1:size
    for y=1:1:size
        A(x,y).left.neighbor2=[[IndexNext(x,size,-1),IndexNext(y,size,1)];[IndexNext(x,size,-1),y];[x,IndexNext(y,size,-1)];[IndexNext(x,size,1),IndexNext(y,size,-1)];[IndexNext(x,size,1),y];[x,IndexNext(y,size,1)]];
        A(x,y).right.neighbor2=[[IndexNext(x,size,-1),IndexNext(y,size,1)];[x,IndexNext(y,size,1)];[IndexNext(x,size,1),y];[IndexNext(x,size,1),IndexNext(y,size,-1)];[x,IndexNext(y,size,-1)];[IndexNext(x,size,-1),y]];
    end
end


for x=1:1:size
    for y=1:1:size
        A(x,y).left.neighbor3=[[IndexNext(x,size,-1),IndexNext(y,size,-1)];[IndexNext(x,size,1),IndexNext(y,size,-1)];[IndexNext(x,size,-1),IndexNext(y,size,1)]];
        A(x,y).right.neighbor3=[[IndexNext(x,size,1),IndexNext(y,size,1)];[IndexNext(x,size,1),IndexNext(y,size,-1)];[IndexNext(x,size,-1),IndexNext(y,size,1)]];
    end
end

end

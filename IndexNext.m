function [next] = IndexNext(now,size,step)


if step>0
    if now+step>size
        next=IndexNext(now-size,size,step);
    else
        next=now+step;
    end
else
    if now+step<1
        next=IndexNext(now+size,size,step);
    else
        next=now+step;
    end
end
end


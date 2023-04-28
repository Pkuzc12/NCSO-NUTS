function [x_next,p_next] = Leapfrog(x,p,epsilon,size,ES,parameter)

p_next=p-epsilon/2*GradEb(x,size,ES,parameter);
x_next=x+epsilon*p_next;
p_next=p_next-epsilon/2*GradEb(x_next,size,ES,parameter);

end


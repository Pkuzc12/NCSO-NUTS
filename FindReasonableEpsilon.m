function [epsilon] = FindReasonableEpsilon(x,kT,size,ES,parameter)

epsilon=1;
p=randn(size*size*4,1);

[x_temp,p_temp]=Leapfrog(x,p,epsilon,size,ES,parameter);

a=2*(Boltzmann(x_temp,kT,size,ES,parameter)/Boltzmann(x,kT,size,ES,parameter)>0.5)-1;

while (Boltzmann(x_temp,kT,size,ES,parameter)/Boltzmann(x,kT,size,ES,parameter))^a>2^(-a)

    epsilon=2^a*epsilon;
    [x_temp,p_temp]=Leapfrog(x,p,epsilon,size,ES,parameter);

end

end


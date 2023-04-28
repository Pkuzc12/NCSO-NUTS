function [Boltzmann] = Boltzmann(x,kT,size,ES,parameter)

Boltzmann=exp(-CalculateEb(x,size,ES,parameter)/(kT));

end


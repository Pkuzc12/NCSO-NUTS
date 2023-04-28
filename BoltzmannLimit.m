function [Boltzmann] = BoltzmannLimit(x,kT,size,ES,parameter,limit)

Boltzmann=exp(limit-CalculateEb(x,size,ES,parameter)/(kT));

end


function [variables] = deriv_vdp(t,vars,mu)
variables=[vars(2)*mu; (1-vars(1)^2)*vars(2)*mu-(vars(1)/mu)];
end


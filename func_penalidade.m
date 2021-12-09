function y=func_penalidade(x,penalidade)
    y = (x(1) - 2)^4 + (x(1) - 2*x(2))^2 + penalidade*(-1/(x(1)^2 - x(2)));
end
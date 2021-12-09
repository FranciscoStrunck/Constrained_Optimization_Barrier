function A=Steepest_Descent(x0, penalidade)
    syms x1 x2
    tol = 1E-6;
    count = 1;
    vet_x{count} = x0;

    while count < 1E3

        g = gradient(func_penalidade([x1, x2], penalidade));
        g = subs(g, [x1; x2], vet_x{count});
        g = double(g);
        d = -g;

        count = count + 1; 
        vet_x{count} = Secaoaurea(d , vet_x{count-1}, penalidade);

        if abs(func_penalidade(vet_x{count},penalidade) - func_penalidade(vet_x{count-1},penalidade)) < tol
            A = vet_x{count};
            break        
        end
    end
end



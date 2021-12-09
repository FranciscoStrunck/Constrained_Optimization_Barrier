function A=Univariante(x0,penalidade)
    tol = 1E-8;
    count = 1;
    vet_x{count} = x0;

    while count < 1E5
        count = count + 1;
        d = [1;0];
        vet_x{count} = Secaoaurea(d , vet_x{count-1},penalidade);

        count = count + 1;
        d = [0;1];
        vet_x{count} = Secaoaurea(d , vet_x{count-1},penalidade);

        if abs(func_penalidade(vet_x{count},penalidade) - func_penalidade(vet_x{count-1},penalidade)) < tol
            A = vet_x{count};
            break
        end
    end
    A = vet_x{count};
end



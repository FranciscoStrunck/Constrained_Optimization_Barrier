function A=Newton(x0, penalidade)
    syms x1 x2
    tol = 1E-8;
    count = 1;
    vet_x{count} = x0;
    
    while count < 1E8      
        
        H = hessian(func_penalidade([x1, x2],penalidade));
        H = subs(H, [x1; x2], vet_x{count});
        S = inv(H);         

        g = gradient(func_penalidade([x1, x2],penalidade));
        g = subs(g, [x1; x2], vet_x{count});

        d = -S*g;
        d = double(d);
        
        count = count + 1; 
        vet_x{count} = Secaoaurea(d , vet_x{count-1}, penalidade);

        B = func_penalidade(vet_x{count},penalidade) - func_penalidade(vet_x{count-1},penalidade);
        if abs(B) < tol
            A = vet_x{count};
            break
        end        
    end    
end







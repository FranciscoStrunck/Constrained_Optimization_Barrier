function A=BFGS(x0, penalidade)
    syms x1 x2
    tol = 1E-4;
    count = 1;
    vet_x{count} = x0;

    % Inicialização e primeira iteração (k=0)

    H = hessian(func_penalidade([x1, x2],penalidade));
    H = subs(H, [x1; x2], vet_x{count});
    H = double(H);
    S = inv(H);
    
    g0 = gradient(func_penalidade([x1, x2],penalidade));
    g0 = subs(g0, [x1; x2], vet_x{count});
    g0 = double(g0);
    g0_transposta = g0.';

    d0 = -S*g0;
    d0_transposta = d0.';

    alfa_num = g0_transposta*d0;
    alfa_den = d0_transposta*H*d0;
    alfa = -alfa_num/alfa_den;

    count = count + 1;

    x = x0 + alfa*d0;
    vet_x{count} = x;

    while count < 1E5
        
        g = gradient(func_penalidade([x1, x2],penalidade));
        g = subs(g, [x1; x2], vet_x{count});
        g = double(g);
        g_transposta = g.';

        delta_x_0 = vet_x{count} - vet_x{count-1};
        delta_x_0_transposta = delta_x_0.';
        
        ga = gradient(func_penalidade([x1, x2],penalidade));
        ga = subs(ga, [x1; x2], vet_x{count});
        ga = double(ga);
        gb = gradient(func_penalidade([x1, x2],penalidade));
        gb = subs(gb, [x1; x2], vet_x{count-1});  
        gb = double(gb);
        delta_g_0 = ga - gb;        

        v1 = gb*(1+alfa*sqrt(abs((d0_transposta*delta_g_0)/(delta_x_0_transposta*gb)))) - ga;
        v1_transposta = v1.';

        w1 = delta_x_0/(delta_x_0_transposta*delta_g_0);
        w1_transposta = w1.';

        
        H = hessian(func_penalidade([x1, x2],penalidade));
        H = subs(H, [x1; x2], x0);
        H = double(H);
        S = inv(H);

        d0 = -(eye(2) + w1*v1_transposta)*S*(eye(2) + v1*w1_transposta)*ga;
        d0_transposta = d0.';

        alfa_num = (ga')*d0;
        alfa_den = d0_transposta*H*d0;
        alfa = -alfa_num/alfa_den;

        count = count + 1;
        d0 = double(d0);
        alfa = double(alfa);
        vet_x{count} =  vet_x{count-1} + alfa*d0;

        B = func_penalidade(vet_x{count},penalidade) - func_penalidade(vet_x{count-1},penalidade);
        if abs(B) < tol
            A = vet_x{count};
            break
        end
    end
end


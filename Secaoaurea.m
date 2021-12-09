function alfaFIM = Secaoaurea(d , x_zero, penalidade)
    e = 1E-6;
    tol = 1E-4;
    iterMAX=1E4;
    D_alfa = 0.01;
    alfa = 1;
    Ra = (sqrt(5)-1)/2;
    count = 1;
    aux=true;

    vetor_x{count}=x_zero;
    GuardaFunc = func_penalidade(x_zero, penalidade);

    if func_penalidade(x_zero-e*d, penalidade) < func_penalidade(x_zero+e*d, penalidade)
        d=-d;
    end

    while aux && count < iterMAX
        count = count+1;
        f1=GuardaFunc;
        
        x=x_zero+alfa*D_alfa*d;

        vetor_x{count}=x;

        GuardaFunc = func_penalidade(x, penalidade);
        f2=GuardaFunc;

        % O novo f(xk+1) é maior que o calculado anteriormente f(xk)?
        if f2 > f1 && f1 ~= f2
            % Se sim, a função NÃO está minimizando mais. Parar o looping.
            aux = false;
        else
            % Se não, a função ESTÁ minimizando. countinuar fazendo um novo passo para o alfa
            alfa=alfa+1;
        end
    end

    % SEÇÃO ÁUREA
    % Intervalor Inicial -> onde está o mínimo
    alfaL = vetor_x{count-1};
    alfaU = vetor_x{count};

    count2 = 1;

    while count2 < iterMAX

        % Tamanho do intervalo
        beta = abs(alfaU-alfaL);

        % Novos pontos
        alfaE = alfaL + (1-Ra)*beta;
        alfaD = alfaL + Ra*beta;

        % Calcular função nos novos pontos
        f1 = func_penalidade(alfaE, penalidade);
        f2 = func_penalidade(alfaE, penalidade);

        if f1 > f2
            % Se f(alfaE) é MAIOR que f(alfaE) -> Quer dizer que o mínimo está
            % entre [alfaL | alfaE]
            alfaL = alfaE;
        else
            % Se f(alfaE) é MENOR que f(alfaE) -> Quer dizer que o mínimo está
            % entre [alfaD | alfaU]
            alfaU = alfaD;
        end

        % Quando alfaU e alfaL são suficientemente próximos -> FIM
        % alfaM é o mínimo
        if abs(alfaU - alfaL) <= tol
            alfaM = (alfaL + alfaU)/2;
            alfaFIM = alfaM;
            break
        end
        count2 = count2+1;
    end
end


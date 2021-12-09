clear all
clc

tic
count = 1;
y_resticao = 0;
k = 0;
rp(count)=10;
beta=0.1;
tol = 1E-8;
penalidade = 10;
penalidade = rp(count);

vet_x{count} = [0;1];

while count < 1E8
     
    % Minimizando   
    count = count + 1;    
    %vet_x{count} = Univariante(vet_x{count-1}, penalidade);
    %vet_x{count} = Powell(vet_x{count-1}, penalidade);
    %vet_x{count} = Newton(vet_x{count-1}, penalidade);
    %vet_x{count} = Steepest_Descent(vet_x{count-1}, penalidade);
    vet_x{count} = Fletcher(vet_x{count-1}, penalidade);
    %vet_x{count} = BFGS(vet_x{count-1}, penalidade);
    %vet_x{count} = Newton_Raphson(vet_x{count-1}, penalidade);
    
    disp('-----------')

    disp('Vetor que Minimiza')
    disp(vet_x{count})
    disp('Valor da Função')
    disp(func_x(vet_x{count}))

    disp('-----------')
    
    % Testar se minimizou a restrição
    if (1/2)*rp(count-1)*(func_restricao(vet_x{count})) < tol
        disp('----------')
        disp('FIM')
        disp('----------')
        break
    end 
    rp(count)=rp(count-1)*beta;   
end

disp('-----------')

disp('Função respeita a restrição?')
if func_restricao(vet_x{count}) <= tol
    disp('Sim')
else
    disp('Não')
end

disp('-----------')

disp('Vetor que Minimiza')
disp(vet_x{count})
disp('Valor da Função')
disp(func_x(vet_x{count}))

disp('-----------')


toc


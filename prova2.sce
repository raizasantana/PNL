clc()
clear()

//Variaveis do problema
alpha = 0;
m = 2;
n = 3;


A=rand(3,1)*rand(1,3);
//A = [1 10; 9 77;8 5]
b = rand(1,3);
b = b';

Z = kernel(A);

xk = [0;0;0];
xk1 = [0;110;0];
lambdak = 0;

//Definicao da f(x)
function fc = f(A,b,x)
    m = (A*x - b)
    fc = m'*m
endfunction

//Gradiente de f(x)
function gr = g(A,b,x)
    gr = 2*(A')*(A*x -b);
endfunction

//Metodo para restrições lineares
//while(Z'*g(A,b,xk1) <> 0)
while(abs(xk1 - xk) < 10^-6)
    disp(xk1);
    disp(xk);
    xk = xk1;
    
    dk = -(Z*Z')*g(A,b,xk);
    
    //Busca linear
    lambda = 1;
    
    while (f(A,b,(xk + lambda*dk)) < (f(A,b,xk) + (alpha*lambda*g(A,b,xk)')*dk)) 
        lambda = 0.1*lambda;
    end
    
    lamdak = lambda;
    xk1 = xk + lambdak*dk;
    
    
end




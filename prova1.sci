clc()
//f(x)
function [fc] = f(A,b,x)
    m = (A*x - b)
    fc = m'*m
endfunction

//Gradiente de f(x)
function gr = g(A,b,x)
    gr = 2*(A')*(A*x -b) ;
endfunction

//Hessiana de f(x)
function hs = h(A)
    hs = 2*(A'*A);
endfunction

//Direcao dfp
function Hk = dfp(alpha,H, A, b, xk)
    xk1 = xk - (alpha*(H*g(A,b,xk)));
    rk = xk1 - xk;
    sk = g(A,b,xk1) - g(A,b,xk);
    Hk = H + (rk*rk')/(rk'*sk) - (H*sk)*(H*sk)'/(sk'*(H*sk))
endfunction

//Direção bfgs
function Hn = bfgs(alpha, H, A, b, xk)
    xk1 = xk - (alpha*(H*g(A,b,xk)));
    rk = xk1 - xk;
    sk = g(A,b,xk1) - g(A,b,xk);
    Hn = (rk*(rk'))*(rk - (H*sk))/(sk'*(rk.^2))

//b = (rk - (H*sk)).(rk*(rk'))/(sk'*(rk.^2))
endfunction

function d = dotVetores(a,b,dim)
    d = 0
    for i = 1:dim
        d = d + a(i,1)*b(i,1)
    end
endfunction


A = [0.78 -0.02 -0.12 -0.14; -0.02 0.86 -0.04 0.06; -0.12 -0.04 0.72 -0.08; -0.14 0.06 -0.08 0.74]
b = [0.76;0.08;1.12;0.68] 
x = [1.6;0.18;0.3;1.2]   


m = 4
n = 4
it = 0

//Variaveis gradiente
xg_temp = [1000;1000;1000;1000]
xkg = x
t = 4
sigma = 0.2;
 alpha = 2;       
 i =0

////Metodo gradiente    
while (sqrt((xkg - xg_temp)'*(xkg - xg_temp)) >= 10^-30)
  
        dirGR = - g(A,b,xkg)   
       
        
        f0 = f(A,b,x)        
        der0 = dotVetores(g(A,b,x),dirGR,n)
        xg_temp = xkg       
        xkg = xkg + t*dirGR

        ft = f(A,b,xkg)
        if(ft > f0 + sigma*t*der0) then
            ft_temp = 0
            while (ft > f0 + sigma*t*der0) or (t>0)
               if(t==0) then
                   break
                  end
                  
                t=t/alpha    
                xkg=xkg+t*dirGR
                ft_temp = ft
                ft=f(A,b,xkg)
                
               
                i = i +1
            end
        else 
            while (ft < f0 + sigma*t*der0)
                if(t==0) then
                    printf("%g\n",i)
                   break
                  end
               
                t=t*alpha    
                xkg=xkg+t*dirGR
                ft=f(A,b,xkg)
                                i = i +1
            end
            t = t/alpha
        end
    it = it + 1
    xkg = xkg + t*dirGR
    ft = f(A,b,xkg)
end
printf("Gradiente: It %d | Fx %g | X", it, ft)
disp(xkg)
   
    
//Variaveis Newton
xn_temp = [100;100;100;100]
xkn = x
t = 2
sigma = 0.2;
alpha = 2;       
i =0

////Metodo Newton    
while (sqrt((xkg - xg_temp)'*(xkg - xg_temp)) > 10^-30)
    
        dirNt = h(A)\-g(A,b,xkn)      
        
        f0 = f(A,b,x)        
        der0 = dotVetores(g(A,b,x),dirNt,n)
        
        xn_temp = xkn       
        xkn = xkn + t*dirNt

        ft = f(A,b,xkn)
        if(ft > f0 + sigma*t*der0) then
            ft_temp = 0
            while (ft > f0 + sigma*t*der0) or (t>0)
               if(t==0) then
                   break
                  end
                  
                t=t/alpha    
                xkn=xkn+t*dirNt
                ft_temp = ft
                ft=f(A,b,xkn)
                
               
                i = i +1
            end
        else 
            while (ft < f0 + sigma*t*der0)
                if(t==0) then
                   break
                  end
                t=t*alpha    
                xkg=xkg+t*dirNt
                ft=f(A,b,xkn)
                i = i +1
            end
            t = t/alpha
        end
    
    it = it + 1
    xkn = xkn + t*dirNt
    ft = f(A,b,xkn)

end
printf("\nNewton: It %d | Fx %g | X", it, ft)
disp(xkg)



//Variaveis BFGS
xb_temp = [100;100;100;100]
h_temp= [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0  1]
xkb = x
t = 1
sigma = 0.2;
alpha = 2;       
i =0    

////Metodo Quase Newton BFGS
while (sqrt((xkg - xg_temp)'*(xkg - xg_temp)) > 10^-30)
  
        H = dfp(alpha,h_temp, A, b, xkb)
        dirBFGS = -H*g(A,b,xkb)      
        h_temp = H
        f0 = f(A,b,x)        
        der0 = dotVetores(g(A,b,x),dirBFGS,n)
        
        xb_temp = xkb
        xkb = xkb + t*dirBFGS
        ft = f(A,b,xkb)
        if(ft > f0 + sigma*t*der0) then
            ft_temp = 0
            while (ft > f0 + sigma*t*der0) or (t >0) //t>0
               if(t==0) then
                   break
                  end
                t=t/alpha    
                xkb=xkb+t*dirBFGS
                ft_temp = ft
                ft=f(A,b,xkb)
                i = i +1
            end
        else 
            while (ft < f0 + sigma*t*der0)
                if(t==0) then
                   break
                  end
                t=t*alpha    
                xkb=xkb+t*dirBFGS
                ft=f(A,b,xkb)
                i = i +1
            end
            t = t/alpha
        end
    
    it = it + 1
    xkb = xkb + t*dirBFGS
    ft = f(A,b,xkb)
end
printf("\nQuase-Newton (BFGS): It %d | Fx %g | X", it, ft)
disp(xkb)


//Variaveis DFP
xd_temp = [100;100;100;100]
h_temp= [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0  1]
xkd = x
t = 1
sigma = 0.2;
alpha = 2;       
i =0    
////Metodo Quase Newton DFP
while (sqrt((xkg - xg_temp)'*(xkg - xg_temp)) > 10^-12)
    if(g(A,b,xkd) == [0;0]) then
        printf("Ponto inicial é um ponto crítico. ")
       xkd = x*t;
       continue
    else
        H = dfp(alpha,h_temp, A, b, xkd)
        dirDFP = -H*g(A,b,xkd)      
        h_temp = H
        f0 = f(A,b,x)        
        der0 = dotVetores(g(A,b,x),dirDFP,n)
        
        xd_temp = xkd       
        xkd = xkd + t*dirDFP
        ft = f(A,b,xkd)
//         printf("f0 %g  ft %g\n",f0 + sigma*t*der0, ft)
        if(ft > f0 + sigma*t*der0) then
            ft_temp = 0
            while (ft > f0 + sigma*t*der0) or ((ft - ft_temp) > 10^-6)
               if(t==0) then
                   break
                  end
                t=t/alpha    
                xkd=xkd+t*dirDFP
                ft_temp = ft
                ft=f(A,b,xkd)
                i = i +1
            end
        else 
            while (ft < f0 + sigma*t*der0)
                if(t==0) then
                   break
                  end
                t=t*alpha    
                xkd=xkd+t*dirDFP
                ft=f(A,b,xkd)
                i = i +1
            end
            t = t/alpha
        end
    
    it = it + 1
    xkd = xkd + t*dirDFP
    ft = f(A,b,xkd)
end
end
printf("\nQuase-Newton (DFP): It %d | Fx %g | X", it, ft)
disp(xkd)


disp(A\b)






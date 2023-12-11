function [x_k2,RelErr,tg] = HeavyBall(g,mu,lamda,nit,tol, n)

x_k0 = g;
x_k1= x_k0;
x_k2 = x_k0;
beta = 0.01;
alpha = 0.01;

[D DT DTD] = DiffOper(sqrt(n));

tg = tic;

RelErr = zeros(nit,1);
for k = 1:nit
    x_old = x_k2;
    y_k2 = x_k1 + beta*(x_k1 - x_k0);
    x_k2 = y_k2 - alpha * (x_k1-g+ 2*lamda*DTD*x_k1);

    x_k0 = x_k1;
    x_k1 = x_k2;

    RelErr(k) = norm(x_k2-x_old,'fro')/norm(x_k2,'fro');

    if RelErr(k) < tol
       break;
    end

end


tg = toc(tg);

       
end

function [B Bt BtB] = DiffOper(N)
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,1) = [];
D(1,1) = 0;
B = [ kron(speye(N),D) ; kron(D,speye(N)) ];
Bt = B';
BtB = Bt*B;
end
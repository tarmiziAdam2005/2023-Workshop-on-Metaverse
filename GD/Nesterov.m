function [x,RelErr,tg] = Nesterov(g,mu,lamda,nit,tol, n)

x = g;
t = 1;
y = x;

[D DT DTD] = DiffOper(sqrt(n));

tg = tic;

for k = 1:nit
    x_old = x;
    x = y - mu*(y-g+ 2*lamda*DTD*y);
    t_old= t;

    t = (k-1)/(k+2);
    y = x + t*(x-x_old);

    RelErr(k) = norm(x-x_old,'fro')/norm(x,'fro');

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
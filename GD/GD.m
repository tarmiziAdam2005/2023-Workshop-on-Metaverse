function [x,RelErr] = GD(g,mu,lamda,nit,tol, n)

x = g;
[D DT DTD] = DiffOper(sqrt(n));

for k = 1:nit
    x_old = x;
    x = x_old - mu*(x-g+ 2*lamda*DTD*x);

    RelErr(k) = norm(x-x_old,'fro')/norm(x,'fro');

    if RelErr(k) < tol
        break;
    end

end
       
end

function [B Bt BtB] = DiffOper(N)
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,1) = [];
D(1,1) = 0;
B = [ kron(speye(N),D) ; kron(D,speye(N)) ];
Bt = B';
BtB = Bt*B;
end
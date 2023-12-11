
function out = amtv(f,lam,beta,tol, Nit)

f = f(:);

u = f;
N         = length(f);


[DD,DDT,DDTDD] = defDDt(N); % 2d Order difference operators.

rho        = zeros(N,1); % Lagrange multiplier

    for k=1:Nit
        
        %Solve the u subproblem
        u = f + DDT*rho;

        %u_old = u;
        %Solve v subproblem via softhresholding
        x    = DD*u - rho/beta;
        v    = shrink(x, lam/beta);

        
        
        %Lagrange multiplier update
        rho   = rho + 0.25*(v - DD*u);
        
        
        %some problem with computing relative error
        %{
        RelErr(k) = norm(u-u_old,'fro')/norm(u,'fro');
        
        if RelErr(k) < tol
            break;
        end
        %} 
    end

    out.sol = u;
    %out.relError = RelErr(1:k);
end


function z = shrink(x,r)
z = sign(x).*max(abs(x)- r,0);
end


function [DD,DDT,DDTDD] = defDDt(N)
%Create a first order difference matrix D
e = ones(N,1);
B = spdiags([e -e], [1 0], N, N);
B(N,1) =  1;

D = B; 
clear B;
DD = D*D;   % Create the 2nd order difference matrix (D^2)
clear D;
% Create the transpose of D
DDT = DD'; %Remember that DT = -D, also called the backward difference.
DDTDD = DD'*DD;
end
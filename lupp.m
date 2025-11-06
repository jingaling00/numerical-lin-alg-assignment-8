% GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING ALGORITHM 6.2
%To solve the n by n linear system for augmented matrix A

for n=1:20
    A=-tril(ones(n,n),-1)+eye(n,n); 
    A(2:n,1)=(n-2:-1:0).'; 
    A(2:n,n)=(n:-1:2).'; 
    A(1,1)=n;A(1,n)=n;
    
    [L,U,p,q] = lup(A);
    A_inv = inv(A);
    L_inv = inv(L);
    U_inv = inv(U);

    A_inf_norm = norm(A,Inf);
    L_inf_norm = norm(L,Inf);
    U_inf_norm = norm(U,Inf);
    
    A_inv_inf_norm = norm(A_inv,Inf);
    L_inv_inf_norm = norm(L_inv,Inf);
    U_inv_inf_norm = norm(U_inv,Inf);
    
    disp('n='), disp(n);
   %{ 
    disp('cond_inf(A)='), disp(A_inf_norm * A_inv_inf_norm);
    disp('cond_inf(L)='), disp(L_inf_norm * L_inv_inf_norm);
    disp('cond_inf(U)='), disp(U_inf_norm * U_inv_inf_norm);
    %}

      
    disp('L='), disp(L);
    disp('U='), disp(U);
    disp('p'), disp(p);
    disp('A^-1='), disp(A_inv);
    disp('L^-1='), disp(L_inv);
    disp('U^-1='), disp(U_inv);
    disp('A_inf_norm'), disp(A_inf_norm);
    disp('L_inf_norm'), disp(L_inf_norm);
    disp('U_inf_norm'), disp(U_inf_norm);           
    disp('A_inv_inf_norm'), disp(A_inv_inf_norm);
    disp('L_inv_inf_norm'), disp(L_inv_inf_norm);
    disp('U_inv_inf_norm'), disp(U_inv_inf_norm); 

end

function [L, U, p] = lup(A)
% LUPP   Gaussian elimination with partial pivoting
%        [L,U,p] = lupp(A)
% returns the LU factorization with permutation vector p
%
% A(p,:) = L * U

N = size(A,1);
L = zeros(N,N);
NROW = 1:N;
NN = N - 1;
ICHG = 0;
OK = true;

for I = 1:NN
    % Step 3: find pivot
    JP = NROW(I);
    AMAX = abs(A(JP,I));
    IMAX = I;
    for IP = I+1:N
        JP = NROW(IP);
        if abs(A(JP,I)) > AMAX
            AMAX = abs(A(JP,I));
            IMAX = IP;
        end
    end
    
    % Step 4: check singularity
    if AMAX <= 1e-20
        OK = false;
        break;
    end
    
    % Step 5: row interchange
    if NROW(I) ~= NROW(IMAX)
        ICHG = ICHG + 1;
        tmp = NROW(I);
        NROW(I) = NROW(IMAX);
        NROW(IMAX) = tmp;
    end
    
    % Step 6–8: elimination
    I1 = NROW(I);
    for J = I+1:N
        J1 = NROW(J);
        XM = A(J1,I)/A(I1,I);
        A(J1,I) = 0;
        A(J1,I+1:end) = A(J1,I+1:end) - XM*A(I1,I+1:end);
        L(J1,I) = XM;
    end
end

% Step 9–11: form L, U, and permutation vector
if OK
    L = L(NROW,:) + eye(N);
    U = A(NROW,:);
    p = NROW;
else
    error('System has no unique solution');
end
end

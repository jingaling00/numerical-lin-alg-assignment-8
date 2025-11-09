for n=1:20
    A=-tril(ones(n,n),-1)+eye(n,n); 
    A(2:n,1)=(n-2:-1:0).'; 
    A(2:n,n)=(n:-1:2).'; 
    A(1,1)=n;A(1,n)=n;
    
    [L,U,p,q] = lupq(A);
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
   
    %disp('cond_inf(A)='), disp(A_inf_norm * A_inv_inf_norm);
    disp('cond_inf(L)='), disp(L_inf_norm * L_inv_inf_norm);
    disp('cond_inf(U)='), disp(U_inf_norm * U_inv_inf_norm);
    

     %{ 
    disp('L='), disp(L);
    disp('U='), disp(U);
    disp('p'), disp(p);
    disp('q'), disp(q);
    disp('A^-1='), disp(A_inv);
    disp('L^-1='), disp(L_inv);
    disp('U^-1='), disp(U_inv);
     %}
    %disp('A_inf_norm'), disp(A_inf_norm);
    
    %disp('L_inf_norm'), disp(L_inf_norm);
    
    %disp('U_inf_norm'), disp(U_inf_norm);           
    %disp('A_inv_inf_norm'), disp(A_inv_inf_norm);
     
    disp('L_inv_inf_norm'), disp(L_inv_inf_norm);
    %disp('U_inv_inf_norm'), disp(U_inv_inf_norm); 
    
end
function [L, U, p, q] = lupq(A)

N = size(A,1);
L = zeros(N,N);
NROW = 1:N;   % Row permutation vector
NCOL = 1:N;   % Column permutation vector
OK = true;

for k = 1:N-1
    % STEP 3: Find the largest pivot element among submatrix
    AMAX = 0;
    imax = k;
    jmax = k;
    for i = k:N
        for j = k:N
            if abs(A(NROW(i), NCOL(j))) > AMAX
                AMAX = abs(A(NROW(i), NCOL(j)));
                imax = i;
                jmax = j;
            end
        end
    end

    % STEP 4: Check for singularity
    if AMAX <= 1.0e-20
        OK = false;
        break;
    end

    % STEP 5: Interchange rows and columns
    if imax ~= k
        temp = NROW(k);
        NROW(k) = NROW(imax);
        NROW(imax) = temp;
    end
    if jmax ~= k
        temp = NCOL(k);
        NCOL(k) = NCOL(jmax);
        NCOL(jmax) = temp;
    end

    % STEP 7-8: Perform Gaussian elimination on permuted matrix
    I1 = NROW(k);
    for i = k+1:N
        J1 = NROW(i);
        XM = A(J1, NCOL(k)) / A(I1, NCOL(k));
        A(J1, NCOL(k)) = 0;
        for j = k+1:N
            A(J1, NCOL(j)) = A(J1, NCOL(j)) - XM * A(I1, NCOL(j));
        end
        L(J1, NCOL(k)) = XM;
    end
end

% STEP 9: Construct L, U, and permutation vectors
if OK
    L = L(NROW, NCOL) + eye(N);
    U = A(NROW, NCOL);
    p = NROW;
    q = NCOL;
else
    error('System has no unique solution (singular matrix).');
end
end

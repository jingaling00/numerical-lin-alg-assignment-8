% GAUSSIAN ELIMINATION WITH BACKWARD SUBSTITUTION ALGORITHM 6.1
% To solve the n by n linear system for augmented matrix A
A = [             7     -35     102     -208     312     -353     292     -158;
         -35     199     -617     1294     -1984     2272     -1898     1034;
         102     -617     1967     -4204     6533     -7554     6358     -3484;
        -208     1294     -4204     9105     -14285     16639     -14088     7756;
         312     -1984     6533     -14285     22577     -26451     22500     -12435;
        -353     2272     -7554     16639     -26451     31135     -26587     14741;
         292     -1898     6358     -14088     22500     -26587     22778     -12664;
        -158     1034     -3484     7756     -12435     14741     -12664     7058];

b =  [  -41;
         265;
        -899;
        2009;
       -3233;
        3842;
       -3309;
        1848];
A_orig = A;
A = [A b];

 single(A) 
 pause 
 N=min(size(A));
 L=zeros(N,N); 
 P=eye(N,N); 
 TRUE = 1;
 FALSE = 0;
 OK=TRUE; 
% STEP 1
% Elimination Process 
 NN = N-1;
 M = N+1;
 ICHG = 0;
 I = 1;
 while OK == TRUE & I <= NN 
% STEP 2
% use IP instead of I
 IP = I;
 while IP <= N & abs(A(IP,I)) <= 1.0e-20 
 IP = IP+1;
 end;
 if IP == M 
 OK = FALSE;
 else
% STEP 3
 if IP ~= I 
 for JJ = 1:M
 C = A(I,JJ);
 A(I,JJ) = A(IP,JJ);
 A(IP,JJ) = C;
 end;
 PI=eye(N,N);
 PI(I,I)=0; PI(IP,IP)=0;
 PI(I,IP)=1; PI(IP,I)=1; 
 P=PI*P; 
 ICHG = ICHG+1;
 ['interchange rows ', num2str(I), ' and ', num2str(IP)]
 single(A)
 pause
 end;
% STEP 4
 JJ = I+1;
 for J = JJ:N
% STEP 5
% use XM in place of l(J,I)
 XM = A(J,I)/A(I,I);
% STEP 6
 for K = JJ:M
 A(J,K) = A(J,K) - XM * A(I,K);
 end;
% Multiplier XM could be saved in A(J,I).
 A(J,I) = 0;
 L(J,I)= XM;
 end;
 end;
 if IP~=I
 L=PI*L; 
 end 
 I = I+1;
 single(A)
 pause
 end;
 if OK == TRUE
% STEP 7
 if abs(A(N,N)) <= 1.0e-20 
 OK = FALSE;
 else
% STEP 8
% start backward substitution
 X=A(:,M)
 pause
 X(N) = A(N,M) / A(N,N);
 X 
 pause
% STEP 9
 for K = 1:NN
 I = NN-K+1;
 JJ = I+1;
 SUM = 0;
 for KK = JJ:N
 SUM = SUM - A(I,KK) * X(KK);
 end;
 X(I) = (A(I,M)+SUM) / A(I,I);
 X 
 pause
 end;
 
 L=L+diag(ones(1,N));
 U=A(:,1:N);
 
 disp('P='), disp(single(P))
 disp('L='), disp(single(L))
 disp('U='), disp(single(U))
 disp('X='), disp(single(X))

 OUP = 1;
 fprintf (OUP, '\n\nwith %d row interchange(s)\n', ICHG);
 end;
 end;
 if OK == FALSE 
 fprintf(1,'System has no unique solution\n');
 end;

 % Addt'l lines for estimating cond number
 r = A_orig * X - b;
 A = [A_orig r];
  single(A) 
 pause 
 N=min(size(A));
 L=zeros(N,N); 
 P=eye(N,N); 
 TRUE = 1;
 FALSE = 0;
 OK=TRUE; 
% STEP 1
% Elimination Process 
 NN = N-1;
 M = N+1;
 ICHG = 0;
 I = 1;
 while OK == TRUE & I <= NN 
% STEP 2
% use IP instead of I
 IP = I;
 while IP <= N & abs(A(IP,I)) <= 1.0e-20 
 IP = IP+1;
 end;
 if IP == M 
 OK = FALSE;
 else
% STEP 3
 if IP ~= I 
 for JJ = 1:M
 C = A(I,JJ);
 A(I,JJ) = A(IP,JJ);
 A(IP,JJ) = C;
 end;
 PI=eye(N,N);
 PI(I,I)=0; PI(IP,IP)=0;
 PI(I,IP)=1; PI(IP,I)=1; 
 P=PI*P; 
 ICHG = ICHG+1;
 ['interchange rows ', num2str(I), ' and ', num2str(IP)]
 single(A)
 pause
 end;
% STEP 4
 JJ = I+1;
 for J = JJ:N
% STEP 5
% use XM in place of l(J,I)
 XM = A(J,I)/A(I,I);
% STEP 6
 for K = JJ:M
 A(J,K) = A(J,K) - XM * A(I,K);
 end;
% Multiplier XM could be saved in A(J,I).
 A(J,I) = 0;
 L(J,I)= XM;
 end;
 end;
 if IP~=I
 L=PI*L; 
 end 
 I = I+1;
 single(A)
 pause
 end;
 if OK == TRUE
% STEP 7
 if abs(A(N,N)) <= 1.0e-20 
 OK = FALSE;
 else
% STEP 8
% start backward substitution
 E=A(:,M)
 pause
 E(N) = A(N,M) / A(N,N);
 E
 pause
% STEP 9
 for K = 1:NN
 I = NN-K+1;
 JJ = I+1;
 SUM = 0;
 for KK = JJ:N
 SUM = SUM - A(I,KK) * X(KK);
 end;
 E(I) = (A(I,M)+SUM) / A(I,I);
 E
 pause
 end;
 
 L=L+diag(ones(1,N));
 U=A(:,1:N);
 
 disp('P='), disp(single(P))
 disp('L='), disp(single(L))
 disp('U='), disp(single(U))
 disp('E='), disp(single(E))

 OUP = 1;
 fprintf (OUP, '\n\nwith %d row interchange(s)\n', ICHG);
 end;
 end;
 if OK == FALSE 
 fprintf(1,'System has no unique solution\n');
 end;

 % Condition number
 cond_approx = (norm(E,Inf)) / (norm(X,Inf) * eps);


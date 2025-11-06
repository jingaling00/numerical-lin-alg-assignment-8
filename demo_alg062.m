% GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING ALGORITHM 6.2
%To solve the n by n linear system for augmented matrix A
A = [ -2.1585e+07    1.0834e+08    -1.6105e+09    3.1580e+09    1.4500e+10    -4.7490e+11    -4.5400e+11    1.9859e+14   ;     
     1.1580e+04    1.2867e+06    6.5520e+06    -9.9350e+07    4.6740e+08    1.0079e+10    -3.7620e+10    -1.2980e+11   ;     
    -3.7530e+02    3.0770e+03    -1.6800e+04    1.5690e+05    5.2670e+06    -9.5700e+06    -6.7120e+08    -1.3853e+10   ;     
    -4.3490e+00    -6.4720e+01    3.8210e+02    1.1500e+02    1.5612e+05    8.1920e+05    2.1256e+07    -9.2120e+07   ;     
     4.4830e-02    5.8390e-01    -6.2560e+00    4.9300e+00    1.0479e+03    -2.6250e+03    8.0280e+04    -3.0930e+05   ;     
     1.7780e-04    -7.8930e-03    -1.1482e-01    1.5853e+00    -5.8850e+00    1.4109e+02    9.3650e+02    2.2140e+03   ;     
    -1.0485e-05    5.4160e-05    6.7940e-04    -4.9540e-03    6.4160e-02    -2.5070e-01    -3.0390e+00    1.1913e+02   ;     
     7.6500e-09    -8.4130e-07    3.9030e-06    -1.0788e-04    4.4820e-04    -1.5652e-02    -1.4059e-01    -5.8270e-01];

b =  [1.97677234255e+14   ;  -1.5696509972e+11   ;  -1.45283601983e+10   ;  -6.9888251969e+07   ;  -2.3059779727e+05   ;  3.2871677648e+03   ;
    1.15900229075e+02   ;  -7.3859861065e-01  ];
A = [A b];

[L,U,p] = lupp(A);

function [L,U,p] = lupp(A)

% single(A)
 N=min(size(A));
 L=zeros(N,N); 
 TRUE = 1;
 FALSE = 0;
 OK=TRUE; 
 M = N+1;
% STEP 1
 for I = 1:N
 NROW(I) = I;
 end;
% initialize row pointer
 NN = N-1;
 ICHG = 0;
 I = 1;
% STEP 2
 while OK == TRUE & I <= NN 
% STEP 3
 JP = NROW(I);
 AMAX = abs(A(JP,I));
 IMAX = I;
 JJ = I+1;
 for IP = JJ:N
 JP = NROW(IP);
 if abs(A(JP,I)) > AMAX 
 AMAX = abs(A(JP,I));
 IMAX = IP;
 end;
 end;
% STEP 4
 if AMAX <= 1.0e-20 
 OK = FALSE;
 else
% STEP 5
% simulate row interchange
 if NROW(I) ~= NROW(IMAX) 
 ICHG = ICHG+1;
 NCOPY = NROW(I);
 NROW(I) = NROW(IMAX);
 NROW(IMAX) = NCOPY;
 %['interchange rows ', num2str(I), ' and ', num2str(IMAX)]
 %single(A(NROW,1:M))
 end;
 I1 = NROW(I);
% STEP 6
 for J = JJ:N
 J1 = NROW(J);
% STEP 7
 XM = A(J1,I)/A(I1,I);
% STEP 8
 for K = JJ:M
 A(J1,K) = A(J1,K)-XM*A(I1,K);
 end;
% Multiplier XM could be saved in A(J1,I)
 A(J1,I) = 0;
 L(J1,I)=XM; 
 %single(A(NROW,1:M))
 end;
 end;
 I = I+1;
 end;
 if OK == TRUE 
  X=A(NROW,M)  ; 
% STEP 9
 N1 = NROW(N);
 if abs(A(N1,N)) <= 1.0e-20 
 OK = FALSE;
% system has no unique solution
 else
% STEP 10
% start backward substitution
 X(N) = A(N1,M) / A(N1,N);
 %X
 % STEP 11
 for K = 1:NN
 I = NN - K + 1;
 JJ = I + 1;
 N2 = NROW(I);
 SUM = A(N2,JJ)*X(JJ);
 for KK = JJ+1:N
 SUM = SUM+A(N2,KK)*X(KK);
 end;
 X(I) = (A(N2,M) - SUM) / A(N2,I);
% X
 end;
 
 L=L(NROW,:)+diag(ones(1,N));
 U=A(NROW,1:N);
 p = NROW;
 %disp('L='), disp(single(L))
 %disp('U='), disp(single(U))
 %disp('X='), disp(single(X))
 
 OUP=1; 
% fprintf (OUP, '\n\nwith %d row interchange(s)\n', ICHG);
% fprintf(OUP, '\nThe rows have been logically re-ordered to:\n');
% for I = 1:N
% fprintf(OUP, ' %2d', NROW(I)); 
% end;
 %fprintf(OUP,'\n');
 %pause
 
 ID=eye(N,N); 
 P=ID(NROW,:);
 % disp('P='), disp(single(P))
 
 end;
 end;
 if OK == FALSE 
 fprintf(1,'System has no unique solution\n');
 end;
end

function [L,U,P,Q] = lucomplete(A) 

% LU razcep matrike A s kompletnim pivotiranjem

% Bor Plestenjak

n = size(A,2);
p = 1:n;
q = 1:n;
for j = 1 : n-1 
   % poiscemo pivotni element
   [maxV, indj] = max(abs(A(j:n,j:n))); 
   [tmp, s] = max(maxV);
   v = indj(s) + j - 1;
   s = s + j - 1;
   % zamenjamo ustrezni vrstici in stolpca
   A([j v],:) = A([v j],:);
   A(:,[j s]) = A(:,[s j]);
   p([j v]) = p([v j]);
   q([j s]) = q([s j]);
   for i = j+1 : n
      A(i,j) = A(i,j)/A(j,j);
      for k = j+1 : n
         A(i,k) = A(i,k) - A(i,j)*A(j,k);
      end
   end
end   
U = triu(A);
L = eye(n)+tril(A,-1);
P = eye(n); P = P(p,:);
Q = eye(n); Q = Q(:,q);

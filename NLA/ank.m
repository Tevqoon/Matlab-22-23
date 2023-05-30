function [l, As] = ank(X, r)
% [l, As] = ank(X, r)
% Po metodi alternirajocih najmanjsih kvadratov najdemo priblizek Y
% ranga r za tenzor X.
% Y najdemo kot vsoto r tenzorjev Y_i ranga 1, ki so definirani s skalarjem l(i)
% ter naborom vektorjev (a_i^1, ... , a_i^d), pri cemer je d red tenzorja X, vektorji
% a_i^j pa so normirani. Shranimo jih v cell arrayu As, za katerega velja
% As{j}(i) = a_i^j.

ns = size(X);
d = length(ns);
As = cell(1, d);
l = zeros(r, 1);
% inicializiramo kot nakljucne matrike
for i = 1:d
    B = rand(ns(i), r);
    [A, l] = normiraj(B);
    As{i} = A;
end
k = d;
koraki = 0;
korakiMax = max(d * 10, 1000);
napakaPrej = 2.0;
napakaZdaj = 1.0;
% napakaZdaj / napakaPrej < 0.99999
while napakaZdaj > 10^-10 && koraki < korakiMax
   napakaPrej = napakaZdaj;
   % dolocimo smer
   if k == d
      k = 1;
   else
      k = k + 1;
   end
   % izracunamo novi Ak
   B0 = resiPrekNormalnegaSistema(X, As, k);
   [A, l] = normiraj(B0);
   As{k} = A;
   % izracunamo napako
   X1 = zeros(size(X));
   for i = 1:r
       palice = cell(1, d);  % i-ti stolpci matrik A1, ... Ad
       for j = 1:d
          palice{j} = As{j}(:, i); 
       end
       X1 = X1 + l(i) * tenzorRanga1(palice);
   end
   % Ne deluje za red > 2: napakaZdaj = norm(X1 - X, 'fro');
   napakaZdaj = (X - X1) .^ 2;
   for i = 1:d
       napakaZdaj = sum(napakaZdaj);
   end
   koraki = koraki + 1;
   fprintf('Koraki: %d Napaka: %f\n', koraki, napakaZdaj);
end
fprintf('Stevilo korakov: %d\n', koraki);

end
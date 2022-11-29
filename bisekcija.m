function [y, korak] = bisekcija(f,a,b,delta, koraki)

% y=BISEKCIJA(f,a,b,delta) bisekcija za funkcijo f na intervalu [a,b], 
% postopek se konca, ko je interval manjsi od delta
%
% Primer uporabe:
% f = @(x) x^3-5*x+1
% bisekcija(f,1,3,1e-10)

% Bor Plestenjak 2020

fa = f(a);  fb = f(b);
e = b-a;
korak = 0;              % nepomembno, stejemo korake za izpis na zaslon
if sign(fa)==sign(fb) 
   disp('Nepravilen interval')
   return
end
while abs(e)>delta
   if korak > koraki
       break
   end
   e = e/2; c = a+e;
   korak = korak+1;           		% nepomembno, stejemo korake za izpis na zaslon
   fc = f(c);
   if sign(fa)==sign(fc)            % dovolj je le preverjanje predznakov
      a = c;  fa = fc;
   else
      b = c;  fb = fc;
   end
end

y = c;
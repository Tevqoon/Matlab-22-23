function [c, r] = regula_falsi(f,a,b,delta, koraki)

fa = f(a);  fb = f(b);
prev = Inf;
c = -Inf; % Inf - (- Inf) = Inf > delta, dokler da pridemo do tretjega koraka
r = 0;      
if sign(fa)==sign(fb) 
   disp('Nepravilen interval')
   return
end
while abs(prev - c) > delta % To bo prvič potencialno True, ko imamo vsaj dve verdnosti
   if r >= koraki % Če dosežemo število korakov, smo končali
       break
   end
   prev = c; % Gledamo razliko med zaporednima približkoma
   c = a + (a - b) * fa / (fb - fa); % abscica abscisnega presečišča sekante
   r = r + 1; % Iteriramo korak 
   fc = f(c);
   if sign(fa)==sign(fc)            % dovolj je le preverjanje predznakov
      a = c;  fa = fc;
   else
      b = c;  fb = fc;
   end
   
end


function y = tgs2(x);

global LL EE;


d = x(1:LL/2);       % Var. der Operatkoeff.
k = x(LL/2+1:end);   % Var. des Terms k*dx/2
c = (1:2:LL-1)';
c2 = c.^2;

% Gleichungen 1 bis L/2 des Systems erstellen
y = [];
for n = 1:LL/2
    % Terme summieren und Gleichungen konstruieren
    y = [y; sum(c2.*d.*sin(c.*k(n)))];
end

% Gleichungen L/2+1 bis L des Systems erstellen
for n = 1:LL/2
    % rechte Seite des Systems: E-Term einbinden
    % Terme summieren und Gleichungen konstruieren
    vorz = -(-1)^(n+LL/2);
    y = [y; sum(c.*d.*cos(c.*k(n)))-1+vorz*EE];
end

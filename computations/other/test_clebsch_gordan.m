%#ok<*AGROW>
 
j = 8;
m = 4;
j1 = 6;
j2 = 6;

m1 = -j1:j1;
m2 = -j2:j2;

n1 = 2 * j1 + 1;
n2 = 2 * j2 + 1;

C1 = sparse(n1 * n2, 1);
C2 = sparse(n1 * n2, 1);

[Cp, m1p, m2p] = clebsch_gordan(j1, j2, j, m);
ind = (m1p + j1) * n2 + m2p + j2 + 1;
C1(ind) = Cp;

disp([Cp, m1p, m2p]);

ind = [];
Cp  = [];
for a = 1:n1
    for b = 1:n2
        C = clebschgordan(j1, m1(a), j2, m2(b), j, m);
        if C ~= 0
            ind = [ind; (m1(a) + j1) * n2 + m2(b) + j2 + 1];
        	Cp  = [Cp; C];
        end
    end
end
C2(ind) = Cp;

disp(['Max difference: ', num2str(max(abs(C1 - C2)))]);

tic;
for a = 1:1000
    Cp = clebsch_gordan(j1, j2, j, m);
end
toc;

tic;
for a = 1:1000
    for b = 1:length(m1p)
        C = clebschgordan(j1, m1p(b), j2, m2p(b), j, m);
    end
end
toc;

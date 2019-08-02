%#ok<*FLUDLR>

% rng(142857);

util_dir = '/Users/srikanthpatala/Dropbox/NCSU_Research/Repos/gb_hsh/gb_hsh_matlab/Util_functions/';
addpath(genpath(util_dir));

w  = 8 * pi * (rand() - 0.5);
th = pi * rand();
ph = 2 * pi * rand();

j = 3;

% Explicit form for j = 1

% ref = [(cos(w / 2) - 1i * sin(w / 2) * cos(th))^2, ...
%        -1i * sqrt(2) * sin(w / 2) * sin(th) * exp(-1i * ph) * (cos(w / 2) - 1i * sin(w / 2) * cos(th)), ...
%        -(sin(w / 2) * sin(th) * exp(-1i * ph))^2;
%        -1i * sqrt(2) * sin(w / 2) * sin(th) * exp(1i * ph) * (cos(w / 2) - 1i * sin(w / 2) * cos(th)), ...
%        1 - 2 * sin(w / 2)^2 * sin(th)^2, ...
%        -1i * sqrt(2) * sin(w / 2) * sin(th) * exp(-1i * ph) * (cos(w / 2) + 1i * sin(w / 2) * cos(th));
%        -(sin(w / 2) * sin(th) * exp(1i * ph))^2, ...
%        -1i * sqrt(2) * sin(w / 2) * sin(th) * exp(1i * ph) * (cos(w / 2) + 1i * sin(w / 2) * cos(th)), ...
%        (cos(w / 2) + 1i * sin(w / 2) * cos(th))^2];
% ref = flipud(fliplr(ref));

U = rotation_mat(j, [w-2*pi, th, ph]);
% U = rotation_mat(j, w-2*pi, th, ph);

R = rotation([cos(ph) * sin(th), sin(ph) * sin(th), cos(th)], w, 2 * j);
R = flipud(fliplr(R));

disp(['Max U R difference: ', num2str(max(abs(U(:) - R(:))))]);

ax = [cos(ph) * sin(th), sin(ph) * sin(th), cos(th)];
q = cos(w/2); Q = ax.*sin(w/2);
a = q-1i*Q(3); lb = -Q(2)-1i*Q(1);

R1 = zeros(2*j+1,2*j+1);
for alp_val = -j:j
    ri = -alp_val + j + 1;
    for al_val = -j:j
        cj = -al_val + j + 1;
        R1(ri, cj) = rotation_wo_svd(a,lb,q, Q, j, alp_val, al_val);
    end
end

disp(['Max U R1 difference: ', num2str(max(abs(U(:) - R1(:))))]);

% tic;
% for a = 1:1000
% %     U = rotation_mat(j, [w, th, ph]);
%     U = rotation_mat(j, w, th, ph);
% end
% toc;
% 
% ax = [cos(ph) * sin(th), sin(ph) * sin(th), cos(th)];
% tic;
% for a = 1:1000
%     R = rotation(ax, w, 2 * j);
% end
% toc;

rmpath(genpath(util_dir));
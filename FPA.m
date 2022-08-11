close all;clear;clc;
tic
% schwefel
fce = 1;

% rastrigin
% fce = 2;

% rosenbrock
% fce = 3;

n_runs = 1000;
p_size = 40;
probswitch = 0.8;
n_iters = 25;
dimensions = 2;


%% bounds, f
if fce == 1
    mesh_step = 10;
    lb = -500 * ones(1,dimensions);
    ub = 500 * ones(1,dimensions);
    [x1,x2] = meshgrid(lb(1):mesh_step:ub(1),lb(1):mesh_step:ub(1));
    f = 2*418.9829 + (-x1.*sin(sqrt(abs(x1))) - x2.*sin(sqrt(abs(x2))));
    fce_name = 'schwefel';
elseif fce == 2
    mesh_step = 0.1;
    lb = -5.12 * ones(1,dimensions);
    ub = 5.12 * ones(1,dimensions);
    [x1,x2] = meshgrid(lb(1):mesh_step:ub(1),lb(1):mesh_step:ub(1));
    f = 20 + x1.^2 + x2.^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2));
    fce_name = 'rastrigin';
elseif fce == 3
    mesh_step = 0.2;
    lb = -5 * ones(1,dimensions);
    ub = 10 * ones(1,dimensions);
    [x1,x2] = meshgrid(lb(1):mesh_step:ub(1),lb(1):mesh_step:ub(1));
    f = 100 * (x2-x1.^2).^2 + (1-x1).^2;
    fce_name = 'rosenbrock';
end

%%
vals = zeros(n_runs,1);best_val = inf;
coords = zeros(n_runs,2);best_coords = [];
seeds = zeros(n_runs,1);best_seeds = zeros(n_runs,1);
cc = 0;

solutions_x = (NaN);
solutions_y = (NaN);
%%
for j = 1:n_runs
    rngseed = round(abs(randi(100000)));
    rng(rngseed);
    s_rng = rng;
    seeds(j,1) = rngseed;

    for i = 1:p_size
      Sol(i,:) = lb + (ub - lb) .* rand(1,dimensions);
      Fitness(i) = Fun(Sol(i,:),fce);
    end

    [fmin,I] = min(Fitness);
    best = Sol(I,:);
    S = Sol; 
    solutions_x(j, 1) = best(1);
    solutions_y(j, 1) = best(2);
    progress_list(j, 1) = Fitness(1);
 
    for t = 1 : n_iters
            for i = 1:p_size                           
               if rand > probswitch
              L = Levy(dimensions);
              dS = L.*(Sol(i,:) - best);
              S(i,:) = Sol(i,:) + dS;
              S(i,:) = simplebounds(S(i,:),lb,ub);
              else
                  epsilon = rand;
                  JK = randperm(p_size);
                  S(i,:) = S(i,:) + epsilon * (Sol(JK(1),:) - Sol(JK(2),:));
                  S(i,:) = simplebounds(S(i,:),lb,ub);
              end

               Fnew = Fun(S(i,:),fce);
                if (Fnew <= Fitness(i))
                    Sol(i,:) = S(i,:);
                    Fitness(i) = Fnew;
               end

              if Fnew <= fmin
                    best = S(i,:);
                    fmin = Fnew;
              end
              solutions_x(j,p_size * (t-1) + i) = best(1);
              solutions_y(j,p_size * (t-1) + i) = best(2);
              progress_list(j,p_size * (t-1) + i) = fmin;
            end
    end
     
     if best_val > fmin
         cc = cc+1;
         best_val = fmin;
         best_coords = best;
         best_seed = rngseed;
         coords(cc,:) = best;
         best_seeds(cc,:) = rngseed;
         vals(cc,:) = fmin;
     end
     progress = [num2str(j) '/' num2str(n_runs) ];
     disp(progress)
end
%%
coords = coords(1:cc,:);
best_seeds = best_seeds(1:cc,:);
vals = vals(1:cc,:);
all_sols_x = [seeds solutions_x];
[xxx, yyy, val_x] = find(all_sols_x == best_seed);
best_run = xxx(1);
best_progress = progress_list(best_run,:);

%% ploty
str_n_eval = ['Pocet evaluaci:' ' ' num2str(n_iters*p_size)];
disp(str_n_eval);
str_best_sol = ['Nejlepsi reseni:' ' ' num2str(best_coords)];
disp(str_best_sol);
str_f_val_min = ['f_val_min =' ' ' num2str(best_val)];
disp(str_f_val_min);

surf(x1, x2, f, 'linestyle',':', 'FaceAlpha', 1);axis tight;hold on 
plot3(best_coords(1), best_coords(2), best_val, 'm.','MarkerSize',20)
plot3(coords(:,1), coords(:,2), vals(:,1), 'k.','MarkerSize',10)
title('Vizualizace nejelpsiho reseni');
hold off

figure;
progress_plot = [best_progress best_progress(end)];
time = (1:(length(progress_plot)));
stairs(time, progress_plot, 'linestyle','-');
title('Vizualizace prubehu nejlepsiho reseni');
toc
%% 
function s = simplebounds(s,lb,ub)
  ns_tmp = s;
  I = ns_tmp < lb;
  ns_tmp(I) = lb(I); 
  J = ns_tmp > ub;
  ns_tmp(J) = ub(J);
  s = ns_tmp;
end

function L=Levy(d)
beta = 3/2;
sigma = (gamma(1+beta) * sin(pi*beta/2) / (gamma((1+beta)/2) * ...
    beta * 2^((beta-1)/2)))^(1/beta);
u = randn(1,d) * sigma;
v = randn(1,d);
step = u./abs(v).^(1/beta);
L = 0.01 * step;
end

function z = Fun(u,fce)
if fce == 1
    % schwefel
    z = 418.9829*length(u) + sum(-u.*sin(sqrt(abs(u))));
elseif fce == 2
    % rastrigin
    z = 10*length(u) + sum(u.^2 - 10*cos(2*pi*u));
elseif fce == 3
    % rosenbrock
    z = (1-u(1))^2+100*(u(2)-u(1)^2)^2;
end
end

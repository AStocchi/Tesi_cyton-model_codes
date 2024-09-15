close all
clear all
clc

j=8 

smooth_data = Fun_lettura_dati();
F_N0 = 2.5765;


%carico dataset nuovo fit cyton_data
load('./nuovo fit cyton data/search_v4_results_points.mat')
start_point = matrix_best;

%punto di partenza dei valori di gamma
resto_vect_gamma = [0.9974    0.9651    0.8128    0.5278    0.2693    0.1187    0.0485    0.0190    0.0073    0.0028];
lun_g = length(resto_vect_gamma);


num_s = size(start_point,1);

%vettor-cell che conterrà tutti gli andamenti dei punti
points_score = {};

%best_conf = zeros(1,8+2+lun_g);
%best_conf(8+2+lun_g) = 10000000;

tic
for z = 1:num_s
    disp('--punto:') 
    z
    
    vect_param = start_point(z,1:8);
    resto_vect_gamma = start_point(z,9:9+lun_g); 
    
    
    %lun_g = length(resto_vect_gamma);
    
    %lancio cyton
    [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit_v2(F_N0, resto_vect_gamma(1), vect_param(2-1), vect_param(3-1), vect_param(8-1), vect_param(9-1), vect_param(4-1), vect_param(5-1), vect_param(6-1), vect_param(7-1), resto_vect_gamma(2:end));
    %evaluation 
    score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);

    %matrice degli score:
    m_score = zeros(1,8+2+lun_g);
    %salvataggio score cyton
    m_score(1,:) = [vect_param, resto_vect_gamma, score];

    
    %%%%vettore delle modifiche: F_g0, F_md0, F_sd0, F_mb0, F_sb0, F_mbn, F_sbn, F_mdn, F_sdn,
    v_counter = zeros(8,1);
    
    %for j = 1:8
        v_counter(j) = (rand()*2-1)*vect_param(j)*0.1;
    %end
    
    vect_param = vect_param + v_counter';
    
    [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit_v2(F_N0, resto_vect_gamma(1), vect_param(2-1), vect_param(3-1), vect_param(8-1), vect_param(9-1), vect_param(4-1), vect_param(5-1), vect_param(6-1), vect_param(7-1), resto_vect_gamma(2:end));
    %evaluation 
    score_n = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
    
    
    %salvataggio score cyton
    m_score(2,:) = [vect_param, resto_vect_gamma, score_n];

    points_score{end+1} = m_score;
    
end
toc

%%
vect_pre = [];
vect_post = [];
for i = 1:num_s
    m_score = points_score{i};
    vect_pre(i) = m_score(1,20);
    vect_post(i) = m_score(2,20);
end

plot(log10(vect_pre))
hold on
plot(log10(vect_post))

%%
vect_diff = vect_pre-vect_post;

plot(log10(abs(vect_diff)))

for i = 1:length(vect_diff)
    if isnan(vect_diff(i))
        vect_diff(i) = 0;
    end
end

media_err = mean(vect_diff);
media_err_abs = mean(abs(vect_diff));

media_err_relat = mean(vect_diff./vect_pre);
media_err_relat_abs = mean(abs(vect_diff)./vect_pre);













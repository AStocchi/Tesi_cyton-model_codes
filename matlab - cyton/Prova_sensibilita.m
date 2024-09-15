

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc


smooth_data = Fun_lettura_dati();
F_N0 = 2.5765;


load('./nuovo fit cyton data/starting_point_for_searchv3.mat')
start_point = m_score1;

num_s = size(start_point,1);

points_score = {};
% points_score_pre = [];  %zeros(1,8+2);
% points_score_post = [];

tic
for z = 1:num_s
    disp('--punto:') 
    z
    
    vect_param = start_point(z,1:(7+2));
    
    %lancio cyton
    [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit(F_N0, vect_param(1), vect_param(2), vect_param(3), vect_param(8), vect_param(9), vect_param(4), vect_param(5), vect_param(6), vect_param(7));
    %evaluation 
    score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
    
    %matrice degli score:
    m_score = zeros(1,8+2);
    %salvataggio score cyton
    m_score(1,:) = [vect_param, score];
    
    %points_score_pre(end+1,:) = m_score(1,:);
    

    %vettore delle modifiche: F_g0, F_md0, F_sd0, F_mb0, F_sb0, F_mbn, F_sbn, F_mdn, F_sdn,
    v_counter = zeros(7+2,1);
    
    for j = 1:9
        v_counter(j) = (rand()*2-1)*vect_param(j)*0.1;
    end
    
    vect_param = vect_param + v_counter';
    
    
    [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit(F_N0, vect_param(1), vect_param(2), vect_param(3), vect_param(8), vect_param(9), vect_param(4), vect_param(5), vect_param(6), vect_param(7));
    %evaluation 
    score_n = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
   
    %salvataggio score cyton
    m_score(2,:) = [vect_param, score_n];
    
    
    %points_score_post(end+1,:) = m_score(2,:);
    
    points_score{end+1} = m_score;
    
end
toc

%%
vect_pre = [];
vect_post = [];
for i = 1:num_s
    m_score = points_score{i};
    vect_pre(i) = m_score(1,10);
    vect_post(i) = m_score(2,10);
end

%%
vect_diff = vect_pre-vect_post;

media_err = mean(vect_diff);
media_err_abs = mean(abs(vect_diff));

media_err_relat = mean(vect_diff./vect_pre);
media_err_relat_abs = mean(abs(vect_diff)./vect_pre);

%codice da far runnare subito dopo la fine di Search_fitting_v2[in questo caso usare solo 8] 
%(o v3[in questo caso 8+2]) per visualizzarne i risultati di convergenza degli score

%load('matrix_best_of_each_run_all_v1.mat')
%load('point_score_all_v1.mat')

matrix_best = zeros(1,8+2+10)%+2+11-1);   %queste righe quando non commentate servono per generare una nuova matrice dei risultati (ALTERNATIVA: se ne può caricare una pre-computata con 'load()' )
tmp = points_score{1};
matrix_best(1,:) = tmp(end,:);   %queste righe
score_run = tmp(:,8+2+10)%+2+11-1);
plot(score_run)
hold on
for i = 2:length(points_score)
    tmp = points_score{i};
    matrix_best(end+1,:) = tmp(end,:);   %queste righe
    score_run = tmp(:,8+2+10)%+2+11-1);
    plot(score_run)    
end


%%
% for i = 1:170
%     matrix_best(i,11)=i;
% end
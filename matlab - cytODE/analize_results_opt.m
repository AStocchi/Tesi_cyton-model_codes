close all
clear all
clc

%load run 1 (50 ripetizioni)
%load('set_score_v2_run22-04-24-n1.mat')

%load run 2 (100 ripetizioni)
%load('set_score_v2_run22-04-24-n2.mat')

%load run 3 (100 ripetizioni - [score*20]^2 )
%load('set_score_v2_run24-04-24-n1.mat')

%LOAD RUN -WORST- DATI KAJAL
% load('./fit-dati-kajal-ode/set_score_v3_run26-06-24-n1-worst.mat')
%LOAD RUN -MID- DATI KAJAL
% load('./fit-dati-kajal-ode/set_score_v3_run27-06-24-n1-mid.mat')
%LOAD RUN -BEST- DATI KAJAL
% load('./fit-dati-kajal-ode/set_score_v3_run27-06-24-n1-best.mat')
%LOAD RUN -BEST- DATI KAJAL
load('./fit-dati-kajal-ode/set_score_v3_run28-06-24-n1-all-ped.mat')

%divido i dati di ogni simulazione per tipo: division parameter, ..., death parameter, score
vettore_score = zeros(length(matrix_scores),1);

trial = matrix_scores{1};
n_gen = length(trial(1,:));

vettori_div = zeros(length(matrix_scores),n_gen);
vettori_birth = zeros(length(matrix_scores),n_gen);
vettori_mor = zeros(length(matrix_scores),n_gen);

for i=1:length(matrix_scores)
    model = matrix_scores{i};
    vettore_score(i)= model(4,1);
    vettori_div(i,:) = model(1,:);
    vettori_birth(i,:) = model(2,:);
    vettori_mor(i,:) = model(3,:);
 
end

%plot degli interi set di parametri
figure(1)
plot(vettori_div');

figure(2)
plot(vettori_birth');

figure(3)
plot(vettori_mor');

figure(4)
plot(vettore_score);

%ricerca del miglior fit
min(vettore_score)
ind_min = find(vettore_score==min(vettore_score))

%calcolo dei fit medi
vettore_div_medio = mean(vettori_div)
vettore_mor_medio = mean(vettori_mor)

%%
%plot dei fit medi e del migliore su tutti gli altri
figure(1)
hold on
plot(vettore_div_medio,'LineWidth',5)
plot(vettori_div(ind_min,:),'LineWidth',2)

figure(3)
hold on
plot(vettore_mor_medio,'LineWidth',5)
plot(vettori_mor(ind_min,:),'LineWidth',2)

%confronto tra medio div e medio morte e best div e best morte
figure(5)
hold on
plot(vettore_div_medio,'LineWidth',5)
plot(vettore_mor_medio,'LineWidth',5)
plot(vettori_div(ind_min,:),'LineWidth',2)
plot(vettori_mor(ind_min,:),'LineWidth',2)

%%
%ordino i modelli per score
sorted_score = sort(vettore_score);

figure(4)
hold on
plot(sorted_score, '-*')

%trovo gli indici dei fit migliori (best 20)
indexes_best_sol = find(vettore_score<=sorted_score(20))

%calcolo il modello medio come prima ma solo sui migliori
vettore_div_best_medio = mean(vettori_div(indexes_best_sol,:))
vettore_mor_best_medio = mean(vettori_mor(indexes_best_sol,:))

%plot della media dei migliori sugli altri
figure(1)
hold on
plot(vettore_div_best_medio,'-*','LineWidth',4)

figure(3)
hold on
plot(vettore_mor_best_medio,'-*','LineWidth',4)

figure(5)
hold on
plot(vettore_div_best_medio,'-*','LineWidth',4)
plot(vettore_mor_best_medio,'-*','LineWidth',4)

%%
%analogo ma (best 5)

%trovo gli indici dei fit migliori (best 5)
indexes_best_sol_5 = find(vettore_score<=sorted_score(5))

%calcolo il modello medio come prima ma solo sui migliori
vettore_div_best_medio_5 = mean(vettori_div(indexes_best_sol_5,:))
vettore_mor_best_medio_5 = mean(vettori_mor(indexes_best_sol_5,:))

%plot della media dei migliori sugli altri
figure(1)
hold on
plot(vettore_div_best_medio_5,'-+','LineWidth',3)

figure(3)
hold on
plot(vettore_mor_best_medio_5,'-+','LineWidth',3)

figure(5)
hold on
plot(vettore_div_best_medio_5,'-+','LineWidth',3)
plot(vettore_mor_best_medio_5,'-+','LineWidth',3)

%%ora bisogna valutare lo score del best, media_20best, media_5best, media


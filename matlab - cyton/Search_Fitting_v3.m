close all
clear all
clc

%non si può usare un metodo "analitico"/"diretto" come il gradient descent
%perchè non posso valutare variazioni della funzione ma solo punti singoli
%quindi uso un'euristica tipo hill-climb. / simulated-anneal. per la
%ricerca stocastica di migliori punti locali
%-->magari anche con passo di aggiornamento variabile nel tempo [che si restringe..] (se capiamo che ha senso?)

%%%pseudo-codice di fitting:
%lancio fun_lettura_dati
%definizione dei parametri del modello (starting point dell'hill climbing)
%lancio fun_cyton_per_fit
%valuto il modello
    %aggiorno i parametri 
    %lancio nuovo cyton
    %valuto aggiornamento e vedo se mantenerlo (volendo possiamo mettere mantenimento probabilistico [con l'idea di andare verso un simulated annealing])
    %->ripeti->

%utilizzo dei dati iniziali presi dal paper cyton-2011
smooth_data = Fun_lettura_dati();
F_N0 = 2.5765;


%load dati worst pedigree
% load('./dati kajal/5_worst_dati_campionati_tempi_corretti.mat')
% smooth_data = campioni;
% F_N0 = 5; 


%fisso il seed a 0
%rng(0)  

%carico dataset con i punti da cui far partire la ricerca 

%set ridotto con 7 starting point dai vari gruppi - con aggiunta di media e
%varianza per processo di morte post-gen.1
% load('try_starting_points_modified.mat')
% start_point = m_score111;

%carico dataset_from kajal_worst_5_pedigree
% load('./dati kajal/fitting dati worst/start_point_for_random_search_v3.mat')
% start_point = m_score11;


%carico dataset nuovo fit cyton_data
load('./nuovo fit cyton data/starting_point_for_searchv3.mat')
start_point = m_score1;

num_s = size(start_point,1);

%vettor-cell che conterrà tutti gli andamenti dei punti
points_score = {};

best_conf = zeros(1,8+2);
best_conf(8+2) = 1000;

for z = 1:num_s
    tic
    disp('--punto:') 
    z
    
    vect_param = start_point(z,1:(7+2));
    
    %guardo modello risultante:
    %vect_param = [1	5.783291	0.100315	4.815692	0.061926	4.710280	0.134035	5.380436	0.533706];
    
    %lancio cyton
    [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit(F_N0, vect_param(1), vect_param(2), vect_param(3), vect_param(8), vect_param(9), vect_param(4), vect_param(5), vect_param(6), vect_param(7));
    %evaluation 
    score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
    
    %matrice degli score:
    m_score = zeros(1,8+2);
    %salvataggio score cyton
    m_score(1,:) = [vect_param, score];
    
    %salvo best configuration (al primo giro viene salvato certamente, per come è definito best_conf
    if score <= best_conf(end,8+2)
        best_conf = m_score(1,:);
    end

    %vettore delle modifiche: F_g0, F_md0, F_sd0, F_mb0, F_sb0, F_mbn, F_sbn
    v_counter = zeros(7+2,1);

    max_step = 750;
    fast_start = 50;
    probab_fissa = 1; %valore fisso di probabilità da superare per salvare il cambiamento quando non è migliorativo
    eps = 0.05;
    
    step = 0;
    while step < max_step

        %definisco le modifiche dei parametri da passare al modello
        if step < fast_start
            v_counter = rand(7+2,1)*2-1;
            if step == fast_start-1
                eps = 0.08;
            end
        else
            v_counter = zeros(7+2,1);
            v_counter(mod(step,7+2)+1) = rand()*2-1; %seleziono un parametro per volta da modificare
            eps = eps - 10^(-4);
        end

        v_F = vect_param + eps*v_counter';
        
        %riporto il primo parametro entro i limiti di validità
        if v_F(1) > 1
            v_F(1) = 1;
        end
        %riporto a zero i valori che diventano negativi (attenzione le medie possono essere negative! 
        %-> non è vero, lo possono essere in quanto distribuzione, ma il significato non lo permette!! )
        v_F(v_F<0) = 0.0001;

        
        %lancio cyton con parametri modificati
        [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit(F_N0, v_F(1), v_F(2), v_F(3), v_F(8), v_F(9), v_F(4), v_F(5), v_F(6), v_F(7));

        %evaluation 
        score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
        
        %salvataggio score cyton (in caso migliori la situazione oppure con una certa probabilità)
        differenza = score - m_score(end,8+2);  %valori negativi indicano un miglioramento (abbassamento) dello score
        if differenza >= 0
            trial = rand(); %un numero casuale tra 0 e 1        
            if trial > probab_fissa + differenza/score
                m_score(end+1,:) = [v_F, score];
                vect_param = v_F;
            end
        else
            %salvataggio di cambiamento migliorativo
            m_score(end+1,:) = [v_F, score];
            vect_param = v_F;        
            if score <= best_conf(end,8+2) %aggiornamento del miglior risultato           
                best_conf =  m_score(end,:);
            end
        end
        step = step + 1;

    end
    
    points_score{end+1} = m_score;
    toc
end



step

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

%utilizzo dati kajal_all_pedigree
%load('./dati kajal/dati_campionati_tempi_corretti.mat')

%load dati worst pedigree
% load('./dati kajal/5_worst_dati_campionati_tempi_corretti.mat')
% smooth_data = campioni;
% F_N0 = 5; 

%fisso il seed a 0
%rng(0)  

%carico dataset con i punti da cui far partire la ricerca 
%load('score_fitting_media_1.mat')
%start_point = m_score1(1,:);

%carico dataset dei 200+ starting point ottenuti dalle grid-search
%load('starting_points.mat')
%start_point = m_score11;

%carico dataset RIDOTTO per test_algoritmo / error_function
%load('try_starting_points.mat')
%start_point = m_score111;

%carico dataset_from kajal_worst_5_pedigree
% load('./dati kajal/fitting dati worst 4 gen/start_point_for_random_search_v2.mat')
% start_point = m_score11;

%carico dataset nuovo fit cyton_data
load('./nuovo fit cyton data/starting_point_for_searchv2.mat')
start_point = m_score1;

num_s = size(start_point,1);

%vettor-cell che conterrà tutti gli andamenti dei punti
points_score = {};

best_conf = zeros(1,8);
best_conf(8) = 1000;

for z = 1:num_s
    tic
    disp('--punto:') 
    z
    
    vect_param = start_point(z,1:7);
    
    %lancio cyton
    [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit(F_N0, vect_param(1), vect_param(2), vect_param(3), vect_param(2), vect_param(3), vect_param(4), vect_param(5), vect_param(6), vect_param(7));
    %evaluation 
    score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
    
    %matrice degli score:
    m_score = zeros(1,8);
    %salvataggio score cyton
    m_score(1,:) = [vect_param, score];
    
    %salvo best configuration (al primo giro viene salvato certamente, per come è definito best_conf
    if score <= best_conf(end,8)
        best_conf = m_score(1,:);
    end

    %vettore delle modifiche: F_g0, F_md0, F_sd0, F_mb0, F_sb0, F_mbn, F_sbn
    v_counter = zeros(7,1);

    max_step = 750;
    fast_start = 50;
    probab_fissa = 1; %valore fisso di probabilità da superare per salvare il cambiamento quando non è migliorativo
    eps = 0.05;
    
    step = 0;
    while step < max_step

        %definisco le modifiche dei parametri da passare al modello
        if step < fast_start
            v_counter = rand(7,1)*2-1;
            if step == fast_start-1
                eps = 0.08;
            end
        else
            v_counter = zeros(7,1);
            v_counter(mod(step,7)+1) = rand()*2-1; %seleziono un parametro per volta da modificare
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
        [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit(F_N0, v_F(1), v_F(2), v_F(3), v_F(2), v_F(3), v_F(4), v_F(5), v_F(6), v_F(7));

        %evaluation 
        score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
        
        %salvataggio score cyton (in caso migliori la situazione oppure con una certa probabilità)
        differenza = score - m_score(end,8);  %valori negativi indicano un miglioramento (abbassamento) dello score
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
            if score <= best_conf(end,8) %aggiornamento del miglior risultato           
                best_conf =  m_score(end,:);
            end
        end
        step = step + 1;

    end
    
    points_score{end+1} = m_score;
    toc
end



step


%% COSE DA FARE:

%->raccogliere i punti salvati in un unico dataset e ciclare su quelli ->
%-> m_score ha ora una dimensione in più (i punti di start)
    %per ogni punto si può far ripetere più volte la run (in particolare i 50 speed_start iniziali) [VEDIAMO SE FARLO DOPO]
        %e salvo solo la run che porta al risultato migliore

%-> speed_start (5-8-10%) delle iteraz. con eps = 0.05 
    %gamma ha dominio [0 1] , %media [1 10], %sigma [0 2.5] -> si può adattare eps di conseguenza se si vuole [VEDIAMO SE FARLO DOPO]

%-> poi dopo speed_start eps=0.075 e lo si porta verso 0.01 (linearmente direi)

    %come aggiunta(o alternativa) [provare su qualche punto CON SEED FISSATO!!] (anche questa solo dopo speed_start??) [VEDIAMO SE FARLO DOPO]
        %far variare la probabilità tipo: 1/2 -> 1/8 ---> 1/4 -> 1/16 ---> 1/8 -> 1/32 -> 1/64 ==> 0.
        %quindi andamento "a scala"/"zig-zag"

    %provare diversi valori fast_start e probab_fissa!! oltre che max_step
        %ovviamente con seed fissato!
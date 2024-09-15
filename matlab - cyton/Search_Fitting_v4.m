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
%load('try_starting_points_modified.mat')

%set di punti alla fine dell'ottimizzazione!
%load('try_optimized_points_modified.mat') 
%start_point = matrix_best;

%set di punti alla fine dell'ottimizzazione v2 o v3 (kajal worst 5 pedigree)
% load('./dati kajal/nuova tranche fitting v3 corretta/worst 4 gen/100_points_from_v2_v3.mat')
% start_point = matrix_best1;

%carico dataset nuovo fit cyton_data
load('./nuovo fit cyton data/result_for_start_searchv4.mat')
start_point = matrix_best1;

%punto di partenza dei valori di gamma
resto_vect_gamma = [0.9974    0.9651    0.8128    0.5278    0.2693    0.1187    0.0485    0.0190    0.0073    0.0028];
lun_g = length(resto_vect_gamma);


num_s = size(start_point,1);

%vettor-cell che conterrà tutti gli andamenti dei punti
points_score = {};

best_conf = zeros(1,8+2+lun_g);
best_conf(8+2+lun_g) = 10000000;

for z = 1:num_s
    tic
    disp('--punto:') 
    z
    
    
    vect_param = start_point(z,1:(7+2));
    
    %prova con risultato migliore ottenuto 
    %vect_param = [0.339459805103434,6.06878202062149,0.860658321970287,6.53803345186243,0.130111530681354,3.77977337760718,0.555196677174966,6.22416185572893,0.484508609393261];
    %risultati ottimizzazione gamma-free +parametri ottimizzati in precedenza (4)
    %vect_param = [0.447899041214090, 6.06878202062149,0.860658321970287,6.53803345186243,0.130111530681354,3.77977337760718,0.555196677174966,6.22416185572893,0.484508609393261]
    %risultati ottimizzazione gamma-free +parametri ottimizzati in precedenza (5)
    %vect_param = [0.806189827170408, 6.23888982688777,1.57528995035417,6.66853960184389,0.341013826039910,4.30796055836781,0.00816646451155349,5.51303121106541,0.000100000000000000]
    
    resto_vect_gamma = [vect_param(1) 0.9974    0.9651    0.8128    0.5278    0.2693    0.1187    0.0485    0.0190    0.0073    0.0028];

%     %prova con risultato migliore ottenuto
%     %resto_vect_gamma = [1,0.943044933527297,0.731010455363113,0.606818980408641,0.429061518432729,0.172237965591217,0.000100000000000000,0.0478235505719339,0.0799990807043874,0.0182844458051435];
    %risultati ottimizzazione gamma-free +parametri ottimizzati in precedenza (4)
    %resto_vect_gamma = [0.447899041214090,0.927504859437558,0.817905177406826,0.706665352774159,0.599269487578909,0.430079661049505,0.175949096308598,0.000100000000000000,0.0131336994925782,0.000100000000000000,0.178657191890817]
    %risultati ottimizzazione gamma-free +parametri ottimizzati in precedenza (5)
    %resto_vect_gamma = [0.806189827170408,0.720210537073639,0.802306797656102,0.732493204115488,0.633090963081080,0.459707437485464,0.193458487884262,0.0376918675173598,0.0981859828833530,0.0757927425127964,0.0228470079688821]
    
    %%%%%%%%%%eeeeeeeeeeeeedit del valore (per variabilità/stabilità del modello)
    %vect_param(2) = vect_param(2)*1.1
    %resto_vect_gamma(8) = resto_vect_gamma(8) - 0.05
    
    
    
    %%%%%%%%%%%%%%%%%%%% valuto punti trovati per all_pedigree
    
    %point 1 -new search tagliata
%     vect_param = [0, 5.88823145542333,0.0441389555349489,4.90967181905796,0.167524080471525,4.66982123234445,0.000100000000000000,5.52689551841923,0.679197232134541]
%     resto_vect_gamma = [1,1,1,0.719670570808332,0.580423601527528,0.225846859869495,0.0434212023577195,0.000100000000000000,0.105382670851385,0.000100000000000000,0.109119561497915,6274605.59625306]
      
    %point 2 -new search tagliata
%     vect_param =[0, 6.29574157551857,0.000100000000000000,4.93570003641201,0.166262477003292,4.65976369988117,0.00336033805967008,5.51927691426852,0.686111963444990]
%     resto_vect_gamma =[1,1,1,0.726695096552460,0.579485702151396,0.353499689935792,0.0956649123360701,0.0302399355324191,0.000100000000000000,0.000100000000000000,0.0230681720406591,6367819.48691950]
     
    
    %point 1
%     vect_param = [0,5.62717725164087,0.632080697186153,4.94421002568102,0.138406974701221,4.65502095059032,0.0516643208938554,5.62717725164087,0.632080697186153]
%     resto_vect_gamma = [1,1,1,0.726870060355567,0.603408788093117,0.205427547803310,0.179376823256418,0.0125476099524422,0.0564769826312214,0.0716247984795582,0.0373921488844041] %7'291'689.980
    
    %point 2
%     vect_param =[0, 6.34621032620005,0.594156475650863,4.74507834126217,0.119197961776262,4.73709599588745,0.0870251670540083,5.76587863102248,0.776757238627836]
%     resto_vect_gamma =[1,1,0.972199706782111,0.728229586523829,0.671223496335953,0.295089553334983,0.168400993231904,0.000938674206122470,0.0608906639942652,0.0472042094331189,0.0660009687584473]  %7'401'903.398
    
    %point 3
%     vect_param =[0, 7.64573011091363,1.51571420014598,4.85545081684303,0.0105402910652787,4.69852585766406,0.138874028620862,5.69521698640427,0.785760990823776]
%     resto_vect_gamma =[1,1,1,0.772080918275313,0.662471461098301,0.156033578637844,0.114672612379351,0.000100000000000000,0.135754443779693,0.191387991164456,0.000100000000000000]  %7'419'188.562
    

        %%%%%%%%%%%%%%%%%%%% valuto best, mid, worst pedigree

        %point 1 -best
%         vect_param = [0,5.38378692803233,0.0786026098888463,4.85170192791364,0.119482145830828,4.66176702916652,0.00375796891718181,5.18970713660982,0.293621269864567]
%         resto_vect_gamma = [1,1,1,0.804462511030055,0.759457303126792,0.0993091086209295,0.205163761243859,0.0590586421554999,0.000100000000000000,0.112545772769985,0.0430823726145989]

        %point 2 -best
%         vect_param = [0,5.82482194314026,0.00553919709520148,4.85410048205943,0.119524903586599,4.66208777729677,0.00585118264502365,5.17540193530589,0.283820662787056]
%         resto_vect_gamma = [1,1,1,0.802936319400436,0.760722092138388,0.365252837464577,0.179868017853272,0.00126145574315776,0.133006985678657,0.127160641951146,0.0394339563274605]

        %point 3 -best
%         vect_param = [0,5.19575197300883,0.000100000000000000,4.86676717003637,0.121338393002058,4.65648942809997,0.00498415157381652,5.19183199877964,0.300210400636026]
%         resto_vect_gamma = [1,1,1,0.812378214445108,0.756577311253041,0.179507629283064,0.000100000000000000,0.0795852597012214,0.0896538281812623,0.000100000000000000,0.000100000000000000]

        %point 1 -mid
%         vect_param = [0,5.41631090687055,0.0603527202863045,4.90985811697834,0.119666687518827,4.66227977939199,0.0783949727355755,5.32644681568464,0.413611112325220]
%         resto_vect_gamma = [1,1,1,0.770624749056090,0.642029915049714,0.176398673747605,0.152244965110705,0.0943112393435072,0.00612622549281557,0.000100000000000000,0.130331736104379]

        %point 3 -mid
%         vect_param = [0,5.61505226933610,0.0970943778822733,4.92136810166465,0.117254245077122,4.65677853279579,0.0837309688406591,5.31533162680373,0.404907478473744]
%         resto_vect_gamma = [1,1,1,0.770337789207286,0.640742773925220,0.181427649347083,0.0856574818516485,0.000100000000000000,0.0577473739265912,0.0108841945550273,0.0178436819984283]

        %point 1 -worst
%         vect_param = [0,5.64301482333357,0.000100000000000000,4.66617165617553,0.000100000000000000,5.04582613431742,0.272224525986023,5.81923895056417,0.739559272780134]
%         resto_vect_gamma = [1,1,1,0.650916991889357,0.347547684657119,0.226315786901488,0.142765872140868,0.000100000000000000,0.129711239715728,0.000100000000000000,0.0920623912146473]
    
        %point 5 -worst
%         vect_param = [0,6.01793767878729,0.109162132535578,4.65161412060009,0.00146563171600711,5.05002157284639,0.273630093742815,5.81866401262654,0.733140673356576]
%         resto_vect_gamma = [1,1,1,0.653795634714233,0.355173679183581,0.187881764857722,0.267149933880125,0.0382562147854449,0.150238992261400,0.0582022815682441,0.0322643453091790]
        
        %point 1 -new worst 4 gen
%         vect_param = [0,5.91908217730045,0.000100000000000000,4.83591475252063,0.333134941880862,4.94252903165055,0.000100000000000000,5.54100126105735,0.288991371674232]
%         resto_vect_gamma = [1,1,0.948207601586732,0.826402171045651,0.494303058611312,0.354059313195992,0.132952323712350,0.0524357828000720,0.0305722759036422,0.0541937477670869,0.00845241193453698]
    
        %point 2 - new worst 4 gen
%         vect_param = [0,5.89104449475015,0.000669060904870625,4.85398318038109,0.325223473172505,4.93412337035299,0.000100000000000000,5.55611245190831,0.296881915262239]
%         resto_vect_gamma = [1,1,0.948479081137582,0.825415255531947,0.561816483957561,0.154929759430114,0.100158078091365,0.0498327819675511,0.0663929191990014,0.000100000000000000,0.0828030443518181]

        %point 5 - new worst 4 gen
%         vect_param = [0,6.00846108685753,0.0450737319200646,4.85916684909092,0.264251491828205,4.92614272928916,0.121511378412339,5.56783841018814,0.281644114262716]
%         resto_vect_gamma = [1,1,0.949212616965383,0.824342538670068,0.490350741648668,0.281219940555577,0.174934423133966,0.0180464383355447,0.0308883880245662,0.0413324826889430,0.0803470803300196]
        %%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%
    
    
    lun_g = length(resto_vect_gamma);
    
    %lancio cyton
    [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit_v2(F_N0, resto_vect_gamma(1), vect_param(2), vect_param(3), vect_param(8), vect_param(9), vect_param(4), vect_param(5), vect_param(6), vect_param(7), resto_vect_gamma(2:end));
    %evaluation 
    score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
    
    %matrice degli score:
    m_score = zeros(1,8+2+lun_g-1);
    %salvataggio score cyton
    m_score(1,:) = [vect_param(2:end), resto_vect_gamma, score];
    
    %salvo best configuration (al primo giro viene salvato certamente, per come è definito best_conf
    if score <= best_conf(end,8+2+lun_g-1)
        best_conf = m_score(1,:);
    end

    %vettore delle modifiche: F_g0, F_md0, F_sd0, F_mb0, F_sb0, F_mbn, F_sbn
    v_counter = zeros(lun_g,1);

    max_step = 1250;
    fast_start = 150;
    probab_fissa = 1; %valore fisso di probabilità da superare per salvare il cambiamento quando non è migliorativo
    eps = 0.05;
    
    step = 0;
    while step < max_step

        %definisco le modifiche dei parametri da passare al modello
        if step < fast_start            
            v_counter = rand(lun_g,1)*2-1;
            if step == fast_start-1
                eps = 0.08;
            end
        else
            v_counter = zeros(lun_g,1);
            v_counter(mod(step,lun_g)+1) = rand()*2-1; %seleziono un parametro per volta da modificare
            eps = eps - 10^(-4);
        end

        
        v_g = resto_vect_gamma + eps*v_counter';
        
        v_F = vect_param;
  
        
        %riporto i parametri entro i limiti di validità (tra 0 e 1) -> visto che sono frazioni.
        v_g(v_g>1) = 1;        
        v_g(v_g<0) = 0.0001;
        

        
        %lancio cyton con parametri modificati
        [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit_v2(F_N0, v_g(1), v_F(2), v_F(3), v_F(8), v_F(9), v_F(4), v_F(5), v_F(6), v_F(7), v_g(2:end));

        %evaluation 
        score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
        
        %salvataggio score cyton (in caso migliori la situazione oppure con una certa probabilità)
        differenza = score - m_score(end,8+2+lun_g-1);  %valori negativi indicano un miglioramento (abbassamento) dello score
        if differenza >= 0
            trial = rand(); %un numero casuale tra 0 e 1        
            if trial > probab_fissa + differenza/score
                m_score(end+1,:) = [v_F(2:end), v_g, score];                
                resto_vect_gamma = v_g;
            end
        else
            %salvataggio di cambiamento migliorativo
            m_score(end+1,:) = [v_F(2:end), v_g, score];
            resto_vect_gamma = v_g;
            if score <= best_conf(end,8+2+lun_g-1) %aggiornamento del miglior risultato           
                best_conf =  m_score(end,:);
            end
        end
        step = step + 1;

    end
    
    points_score{end+1} = m_score;
    toc
end



step

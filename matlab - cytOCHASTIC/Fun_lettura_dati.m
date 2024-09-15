% close all
% clear all
% clc

function [smooth_data] = Fun_lettura_dati()

    raw_data = {};
    data_gen = zeros(9,1);

    %for di lettura dei dati raw dai file .csv
    for i = 0:8
        file_n = "screenshot_cyton_data\Default Dataset(" + i + ").csv";
        data_gen = readmatrix(file_n, "Delimiter",";");
        data_gen(data_gen<0) = 0;
        raw_data{end+1} = data_gen;
        AAflag_plot_raw = 0;
        if AAflag_plot_raw == 1
            plot(data_gen(:,1),data_gen(:,2),'-*');
            tt = "raw gen " + i;
            title(tt)
            hold on    
        end
    end


    %for di smoothing dei dati:
    all_data = {};
    for i = 1:9 %considero un generazione alla volta
        editing = raw_data{i};
        edited = {};
        len = length(editing);
        temp = zeros(1,2);
        for j = 1:len
            if j==1
                %primo punto viene sempre incorporato
                temp(1,:) = editing(j,:);
            else
                diff = abs(editing(j,1)-editing(j-1,1));
                if diff>=1
                    %caso in cui sono passato a diverso istante temp.
                    %faccio media dei punti al "tempo precedente" e poi passo al nuovo punto
                    smo_x = mean(temp(:,1)); %median(temp(:,1)); 
                    smo_y = mean(temp(:,2)); %median(temp(:,2)); 
                    edited{end+1} = [smo_x, smo_y];  %salvo il punto smoothed
                    temp = zeros(1,2);
                    temp(1,:) = editing(j,:);
                else
                    %caso in cui sto guardando lo stesso istante temporale (a meno di errore di misurazione)
                    temp(end+1,:) = editing(j,:);                
                end
                if j == len  %se sono alla fine devo salvare l'istante temporale (anche se non ne trovo uno superiore [che ricade nel caso del if precedente])
                    smo_x = mean(temp(:,1)); %median(temp(:,1));
                    smo_y = mean(temp(:,2)); %median(temp(:,2));
                    edited{end+1} = [smo_x, smo_y];  %salvo il punto smoothed
                end
            end                
        end
        all_data{end+1} = edited;
    end

    %ciclo per trasformare i dati da cells di cell-array ---to--> cells di matrici
    smooth_data = {};
    for i = 1:length(all_data)
        temp = all_data{i};   
        for j = 1:length(temp)
            if j == 1
                te_mp = temp{j};
            else
                te_mp(end+1,:) = temp{j}; 
            end
        end
        smooth_data{end+1} = te_mp;
    end

    AAflag_plot_smoothed = 0;
    if AAflag_plot_smoothed == 1
        for i = 1:length(smooth_data)
            %plot delle singole generazioni smoothed
            gen = smooth_data{i};
            plot(gen(:,1),gen(:,2),'-*');
            tt = "smoothed gen " + (i-1);
            title(tt)
            hold on 
        end
    end


    %ciclo per calcolare la somma totale di cellule ad ogni istante di tempo -> somma delle generazioni
    total_cell = [];
    for i = 1:length(smooth_data)
        gen = smooth_data{i};
        if i == 1
            total_cell = gen;
        else
            %total_cell(:,1) = total_cell(:,1) + gen(:,1);
            total_cell(:,2) = total_cell(:,2) + gen(:,2);
        end
    end
    
    if AAflag_plot_smoothed == 1
        %plot del numero di cellule totali nel tempo
        plot(total_cell(:,1),total_cell(:,2),'-*');
        tt = "smoothed total cell in time";
        title(tt)
    end


    %smooth_data  %dati più puliti che sono riuscito a creare
    
    %% codice di prova per duplicare il numero di punti del segnale
%     augm_data = {}
%     for i = 1:length(smooth_data)
%         tmp = smooth_data{i};
%         augm_signal = zeros(size(tmp,1)*2-1,size(tmp,2));
%         for j = 1:size(augm_signal,1)
%             if mod(j,2)==1
%                 %per j dispari: inserisco il punto originale estratto dai dati
%                 augm_signal(j,:)=tmp((j+1)/2,:);
%             else
%                 %per j pari: faccio la media tra i due punti consecutivi nei dati
%                 augm_signal(j,:)=(tmp((j)/2,:)+tmp((j+2)/2,:))./2;
%             end
%         end
%         augm_data{end+1}=augm_signal;
%     end
%     
%     augm_data
%     
%     smooth_data = augm_data
    
       
end


% close all
% clear all
% clc

function [score] = Fun_evaluation_fit(smooth_data, t_ret, ncell_ret)
        
    score = 0;
    
    num_gen = size(ncell_ret,1);
    
    for i = 1:num_gen
        
        gen = smooth_data{i};
        
        AAflag_plot_complete = 0;
        if AAflag_plot_complete
            %figure(i)
            plot(t_ret(1:end-1),ncell_ret(i,:),'-o')
            hold on            
            plot(gen(:,1),gen(:,2),'-*');
            %hold off
        end
        
        %vettore per salvare gli indici degli istanti di tempo della simulazione da confrontare
        index_t_istant_confronto = zeros(length(gen(:,1)),1);

        %%%ricerca per ogni tempo dei dati del valore della simulazione 
        j = 1;   %contatore sugli istanti dei dati
        t_f = gen(j,1);
        soglia_min = 1;
        m = soglia_min;
        %ricerca LINEARE del tempo della simulazione più vicino da confrontare con i dati 
        k = 1;  %contatore sugli istanti simulativi
        while k<=length(t_ret)
            if abs(t_f-t_ret(k)) < m
                index_t_istant_confronto(j) = k;
                m = abs(t_f-t_ret(k));
            elseif m~=soglia_min
                %puoi entrare qui solo quando la distanza dei tempi è
                %maggiore della soglia (quindi la condizione del I if non è rispettata)
                %ma 'm' deve essere != 10 -> devi essere passato per forza dentro il I if -->
                %--> hai ricominciato ad allontanarti nel tempo (hai superato il più vicino) 

                m = soglia_min;  %setto nuovamente la soglia al massimo
                j = j +1; 
                if j > length(index_t_istant_confronto)
                    break
                end
                t_f = gen(j,1); % e passo a cercare l'istante successivo

                k = k-1;
                %'k=k-1' -> operazione necessaria ad evitare che si salti il controllo del primo passo simulativo dopo uno di quelli scelti 
                %(ci permette di gestire il caso di due istanti consecutivi entrambi da confrontare [sappiamo bene però che è un caso ultra remoto vista la risoluzione temporale scelta])
            end
            k = k+1;
        end
        
        AAflag_plot_confront_signal = 0;
        conf_sign = ncell_ret(i,index_t_istant_confronto);
        if AAflag_plot_confront_signal  
            figure(i)
            %plot del segnale di fitting simulato
            plot(t_ret(index_t_istant_confronto),conf_sign,'-o')
            hold on
            gen = smooth_data{i};
            %plot dei dati
            plot(gen(:,1),gen(:,2),'-*');
            
            
            hold off
        end
        
        %PEZZO DI FUNZIONE CHE SERVE A VALUTARE IL SEGNALE SOLO DOVE IL
        %RIFERIMENTO è >0 (anzi >0.001)
        rif_sign = gen(:,2);
        differenza_finale = rif_sign-conf_sign';
        differenza_finale(rif_sign<0.001)=0;  %azzero il valore dell'errore dove il segnale riferimento è <
        
        %aggiorno lo score aggiungendo il valore assoluto della distanza tra i dati e il segnale simulato
        score = score + sum((abs(differenza_finale).*20).^2);
        
        
                       
    end
    %score   %
end

function [res] = fsigm(vect)
    tmp = vect; 
    tmp = tmp*10;
    tmp = tmp - 8;
    expon = exp(-tmp);
    res = 1 ./ (expon + 1);
end
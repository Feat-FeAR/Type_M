% https://it.mathworks.com/matlabcentral/answers/248797-i-need-to-convert-a-number-into-its-column-name-equivalent

function xlcol_addr=num2xlcol(col_num)

    % Convert numbers to Excel-style column indexes
    % col_num - positive integer greater than zero
    
    n=1;
    while col_num>26*(26^n-1)/25
        n=n+1;
    end
    base_26=zeros(1,n);
    tmp_var=-1+col_num-26*(26^(n-1)-1)/25;
    for k=1:n
        divisor=26^(n-k);
        remainder=mod(tmp_var,divisor);
        base_26(k)=65+(tmp_var-remainder)/divisor;
        tmp_var=remainder;
    end
    xlcol_addr=char(base_26); % Character vector of xlcol address
end

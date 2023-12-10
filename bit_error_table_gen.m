function [bit_error_table] = bit_error_table_gen()
%bit_error_table_gen: 各数字对应的误比特数

M = 4;
bit_error_table = zeros(2^M,2^M);
data = 0:1:2^M-1;
bin_data = dec2bin(data);

for i = 1:1:2^M
    for j=i:1:2^M
        error_num = 0;
        for xx =1:1:M
            if bin_data(i,xx)~= bin_data(j,xx)
                error_num = error_num + 1;
            end
        end
        bit_error_table(i,j) = error_num;
        bit_error_table(j,i) = error_num;
    end
end


end


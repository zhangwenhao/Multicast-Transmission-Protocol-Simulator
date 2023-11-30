function [r_1]=fun_db_to_math(xdb)
      %xdb=10*log10(r_1);
      r_1=power(10,xdb/10);
end
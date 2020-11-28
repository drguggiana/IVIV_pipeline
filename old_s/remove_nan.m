function str_o=remove_nan(str,temp)

emptyIndex = find(arrayfun(@(str) isempty(temp),str))
for i=1:length(emptyIndex)
   temp(empty_Index(i))=NaN;
 
end
str_o=str;
end

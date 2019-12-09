function WriteWorkfile(OriginalFile, oldstring, newstring)

fin = fopen(OriginalFile,'rt');
fout = fopen('WorkFile.f90','wt');

while ~feof(fin)
   s = fgetl(fin);
   s = strrep(s, oldstring, newstring);
   fprintf(fout,'%s\n',s);
   % disp(s)
end

fclose(fin);
fclose(fout);
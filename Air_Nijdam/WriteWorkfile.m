function WriteWorkfile(OriginalFile, Outputfile, oldstring, newstring)

fin = fopen(OriginalFile,'rt');
fout = fopen(Outputfile,'wt');

while ~feof(fin)
   s = fgetl(fin);
   s = strrep(s, oldstring, newstring);
   fprintf(fout,'%s\n',s);
   % disp(s)
end

fclose(fin);
fclose(fout);
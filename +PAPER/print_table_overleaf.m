function [src] = print_table_overleaf(src)
direc = '~/projects/vis/report/';

ffile = [direc 'report.tex'];
fid = fopen(ffile,'wt');
fprintf(fid, src);
fclose(fid);

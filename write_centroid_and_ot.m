function [] = write_centroid_and_ot(c, OT, prefix_name)

fp=fopen([prefix_name, '_c.d2'], 'w');
for i=1:length(c)
    fprintf(fp, '%d\n', size(c{i}.supp,1));
    fprintf(fp, '%d\n', size(c{i}.supp,2));
    fprintf(fp, '%f ', c{i}.w);
    fprintf(fp, '\n');
    for j=1:length(c{i}.w)
        fprintf(fp, '%f ', c{i}.supp(:,j));
        fprintf(fp, '\n');
    end
end
fclose(fp);

fp=fopen([prefix_name, '_c.ot'], 'w');
for i=1:length(c)
    for j=1:length(OT{i})
        fprintf(fp, '%f ', OT{i}{j}(:));
        fprintf(fp, '\n');
    end
end
fclose(fp);

end
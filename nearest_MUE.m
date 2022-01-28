function d = nearest_MUE(xFBS, yFBS, mue)
    d = inf;
for i=1:size(mue,2)
    d_next = sqrt((xFBS-mue(i).X)^2+(yFBS-mue(i).Y)^2);
    if d_next < d
        d = d_next;
    end
end
end
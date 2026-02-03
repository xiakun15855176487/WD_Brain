function genesymbolOv = genesymbolOverlap(genesymbol1,genesymbol2)
genesymbolOv = [];
for j = 1:length(genesymbol1)
    for k = 1:length(genesymbol2)
        if strcmp(genesymbol1(j),genesymbol2(k))
            genesymbolOv = [genesymbolOv;genesymbol1(j)];
            break
        end
    end
end
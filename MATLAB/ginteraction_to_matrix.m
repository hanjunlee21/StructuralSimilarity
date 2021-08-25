clear all;

filenames = ["GM12878_1_2.chr16.25kb.KR.obs_exp";"K562.chr16.25kb.KR.obs_exp"];

parfor filecnt = 1:size(filenames,1)
    tsv = readmatrix(strcat("../AggregateContacts/",filenames(filecnt,1),".tsv"), 'FileType', 'text', 'Delimiter', '\t', 'OutputType', 'double');
    
    regions = [62.5,78];
    
    window = 2e6;
    step = 500e3;
    bin = 25e3;
    
    for regioncnt = 1:size(regions,1)
        
        output = ones(window/bin,window/bin);
        
        % diagonal
        diagonal = tsv((tsv(:,6)-tsv(:,3))==0,:);
        for j = 1:window/bin
            if sum(diagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin) ~= 0
                output(j,j) = diagonal(diagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin,7);
            end
        end
        
        % offdiagonal
        for k = 1:(window/bin-1)
            offdiagonal = tsv((tsv(:,6)-tsv(:,3))==k*bin,:);
            
            for j = 1:window/bin-k
                if sum(offdiagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin) ~= 0
                    output(j,j+k) = offdiagonal(offdiagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin,7);
                    output(j+k,j) = output(j,j+k);
                end
            end
        end
        
        writematrix(output, strcat("../Regions/",filenames(filecnt,1),".",num2str(regions(regioncnt,1)),"Mb.tsv"), 'FileType', 'text', 'Delimiter', '\t');

    end
end

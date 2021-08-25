clear all;

% Region (Mb)
regions = [62.5;78];

filenames = ["GM12878_1_2.chr16.25kb";"K562.chr16.25kb"];
tsv1 = readmatrix(strcat("../AggregateContacts/",filenames(1,1),".tsv"), 'FileType', 'text', 'Delimiter', '\t', 'OutputType', 'double');
tsv2 = readmatrix(strcat("../AggregateContacts/",filenames(2,1),".tsv"), 'FileType', 'text', 'Delimiter', '\t', 'OutputType', 'double');

filenames = ["GM12878_1_2.chr16.25kb.KR.obs_exp";"K562.chr16.25kb.KR.obs_exp"];

for regioncnt = 1:size(regions,1)
    mat1 = readmatrix(strcat("../Regions/",filenames(1,1),".",num2str(regions(regioncnt,1)),"Mb.tsv"), 'FileType', 'text', 'Delimiter', '\t', 'OutputType', 'double');
    mat2 = readmatrix(strcat("../Regions/",filenames(2,1),".",num2str(regions(regioncnt,1)),"Mb.tsv"), 'FileType', 'text', 'Delimiter', '\t', 'OutputType', 'double');
    
    window = 2e6;
    step = 500e3;
    bin = 25e3;
    k1 = 0.01;
    k2 = 0.03;
    L = 1;
    c1 = (k1*L)^2;
    c2 = (k2*L)^2;
    c3 = c2/2;
    ssimwindow = 7;
    
    luminance = ones(window/bin+1-ssimwindow,window/bin+1-ssimwindow);
    contrast = ones(window/bin+1-ssimwindow,window/bin+1-ssimwindow);
    structure = ones(window/bin+1-ssimwindow,window/bin+1-ssimwindow);
    for x = 1:window/bin+1-ssimwindow
        for y = 1:window/bin+1-ssimwindow
            window1 = reshape(mat1(x:x+ssimwindow-1,y:y+ssimwindow-1),ssimwindow*ssimwindow,1);
            window2 = reshape(mat2(x:x+ssimwindow-1,y:y+ssimwindow-1),ssimwindow*ssimwindow,1);
            
            %luminance
            luminance(x,y) = (2*mean(window1)*mean(window2)+c1) / (mean(window1)*mean(window1)+mean(window2)*mean(window2)+c1);
            luminance(y,x) = luminance(x,y);
            
            %contrast
            contrast(x,y) = (2*std(window1)*std(window2)+c2) / (std(window1)*std(window1)+std(window2)*std(window2)+c2);
            contrast(y,x) = contrast(x,y);
            
            %structure
            covariancematrix = cov([window1, window2]);
            covariance = covariancematrix(1,2);
            structure(x,y) = (covariance+c3) / (std(window1)*std(window2)+c3);
            structure(y,x) = structure(x,y);
        end
    end
    ssim = luminance.*contrast.*structure;
    
    [~,minx] = min(min(ssim).');
    [~,minys] = min(ssim);
    miny = minys(1,minx);
    
    disp(strcat(num2str(regions(regioncnt,1))," Mb"));
    disp([minx, miny]);
    
    window1 = reshape(mat1(minx:minx+ssimwindow-1,miny:miny+ssimwindow-1),ssimwindow*ssimwindow,1);
    window2 = reshape(mat2(minx:minx+ssimwindow-1,miny:miny+ssimwindow-1),ssimwindow*ssimwindow,1);
    
    %luminance
    l = (2*mean(window1)*mean(window2)+c1) / (mean(window1)*mean(window1)+mean(window2)*mean(window2)+c1);
    disp(strcat("Luminance: ", num2str(l)));
    
    %contrast
    c = (2*std(window1)*std(window2)+c2) / (std(window1)*std(window1)+std(window2)*std(window2)+c2);
    disp(strcat("Contrast: ", num2str(c)));
    
    %structure
    covariancematrix = cov([window1, window2]);
    covariance = covariancematrix(1,2);
    s = (covariance+c3) / (std(window1)*std(window2)+c3);
    disp(strcat("Structure: ", num2str(s)));
    
    disp(strcat("SSIM: ", num2str(l*c*s)));
    
    % import coverage data
    
    coverage1 = zeros(window/bin,window/bin);
    
    % diagonal
    diagonal = tsv1((tsv1(:,6)-tsv1(:,3))==0,:);
    for j = 1:window/bin
        if sum(diagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin) ~= 0
            coverage1(j,j) = diagonal(diagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin,7);
        end
    end
    
    % offdiagonal
    for k = 1:(window/bin-1)
        offdiagonal = tsv1((tsv1(:,6)-tsv1(:,3))==k*bin,:);
        
        for j = 1:window/bin-k
            if sum(offdiagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin) ~= 0
                coverage1(j,j+k) = offdiagonal(offdiagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin,7);
                coverage1(j+k,j) = coverage1(j,j+k);
            end
        end
    end
    
    coverage2 = zeros(window/bin,window/bin);
    
    % diagonal
    diagonal = tsv2((tsv2(:,6)-tsv2(:,3))==0,:);
    for j = 1:window/bin
        if sum(diagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin) ~= 0
            coverage2(j,j) = diagonal(diagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin,7);
        end
    end
    
    % offdiagonal
    for k = 1:(window/bin-1)
        offdiagonal = tsv2((tsv2(:,6)-tsv2(:,3))==k*bin,:);
        
        for j = 1:window/bin-k
            if sum(offdiagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin) ~= 0
                coverage2(j,j+k) = offdiagonal(offdiagonal(:,3)==regions(regioncnt,1)*1e6+(j-1)*bin,7);
                coverage2(j+k,j) = coverage2(j,j+k);
            end
        end
    end
    
    save(strcat("../Figure2/",num2str(regions(regioncnt,1)),"Mb.mat"),'ssim','window1','window2','coverage1','coverage2');
end


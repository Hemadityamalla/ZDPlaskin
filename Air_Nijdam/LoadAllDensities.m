clear;

filebase = 'D:\Werk\ZDPlasKin\Sander\Output\Pressure\';

%OxygenConc = [0 1e-6 1e-5 1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1 2e-1];
Pressures = [10 25 50 75 100 133 150 180 210 250 300];
edens_initial = 1e9.*power(Pressures./133,1.66);

for i=1:size(Pressures,2)
    % foldernames{i} = sprintf('%g',OxygenConc(i));
    foldernames{i} = sprintf('%g',Pressures(i));
end

GroupTitles = {'All'};

for i=1:size(foldernames,1)
    for j=1:size(foldernames,2)
        filenames{i,j} = [filebase foldernames{i,j} '\qt_densities.txt'];
        [Times{i,j},EDens{i,j}]=importedensities(filenames{i,j});
    end
end

% ProbeDens = [4E8 3E8 2E8 1E8 7E7]*2;
DoPercentages = true;
ProbeDens = [0.01 0.1 0.5 0.8];
PulseLength = 0;

PlotDensityTimes;

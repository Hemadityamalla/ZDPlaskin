OxygenConc = [0 1e-6 1e-5 1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1 2e-1];
%OxygenConc = [3e-3];

Pressures = [10 25 50 75 100 133 150 180 210 250 300];

for k=1:size(Pressures,2);
    
    OutputFolder = sprintf('Output\\Pressure\\%g',Pressures(k));
    
%    WriteNewOxygen('Sander-Laser-Decay-WorkFile.f90',OxygenConc(k));
%    WriteNewOxygen('Aram_N2_O2.f90', 'oxygen_conc = 0.01', sprintf('oxygen_conc = %0.3e',OxygenConc(k)));
    WriteWorkfile('Sander-Laser-Decay-Pressure-WorkFile.f90','Temp.f90',...
        'pressure = 133.0d-3', sprintf('pressure = %3.1dd-3',Pressures(k)));
    
    edens(k) = 1e9*power(Pressures(k)/133,1.66);
    
    WriteWorkfile('Temp.f90','WorkFile.f90',...
        'e_density = 1e9', sprintf('e_density = %0.3e',edens(k)));

    disp(['Working file number ' num2str(k)]);
    system('gfortran dvode_f90_m.F90 zdplaskin_m.F90 Workfile.F90 bolsig_g.dll');
    system('a');
    
    mkdir(OutputFolder);
    copyfile('WorkFile.f90',OutputFolder);
    movefile('qt*',OutputFolder);
end
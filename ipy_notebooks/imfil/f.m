import matlab.io.hdf4.*

format long
 
% This script is designed to bridge between imfil.m and sample_obj_function.py

% Set options if desired, have not seen them have an effect
options=imfil_optset;
% options=imfil_optset('termtol',10,options);
% options=imfil_optset('function_delta',-20,options);
% options=imfil_optset('verbose',1,options);
% options=imfil_optset('scaledepth',10,options);
options=imfil_optset('scalestart',1,options);
% options=imfil_optset('fscale',-1.01,options);
% options=imfil_optset('random_stencil',2,options);
% options=imfil_optset('add_new_directions','tangent_directions',options);


% Read in initial values
ampID = fopen('1.57rad_init_amp.txt');
initamps = [];

% init = strrep(fgetl(ampID),';',' ');
 init = fgetl(ampID);
while ischar(init)
    initvec = reshape(str2num(char(strsplit(init))),2,[]);
    initamps = [initamps initvec];
%     init = strrep(fgetl(ampID),';',' ');
    init = fgetl(ampID);
end

fclose(ampID);

disp(initamps)
outname='imfil_1.57rad_out_0_noise_mod.csv';

initamps2 = [];
results = [];
fvalues = [];
iterations = [];

%Try multiple initial values 
n=1;
while n < 11
    initvec = [initamps(1,n); initamps(2,n)];
    disp(initvec);

    %Inner loop feeds results of previous minimization to check for
    %convergence
    inner=1;
    while inner < 3
    
        initamps2 = [initamps2; initvec.'];

        [result, histout, comphist] = imfil(initvec, @func, 30, [-1. 1.; -1. 1.], options);

        results = [results; result.'];
        fvalues = [fvalues; histout(end, 2)];
        iterations = [iterations; histout(end,1)];


        initvec = result;
        disp(initvec);
        inner = inner+1;
    end;
    
    n = n+1;
end;

disp("--")
disp(initamps2)
disp("--")
disp(results)
disp("--")
disp(fvalues)
disp("--")
disp(iterations)
dlmwrite(outname, [initamps2 results fvalues iterations], 'delimiter', ',', 'precision', 10); 

% fclose(outID);

function [fout, ifail, icount] = func(x,h)
    % Run external python objective function here, through command line
    % integrated into matlab. Include x as cmd line arguments
    
    %!C:\Users\hp\Documents\GitHub\CircuitNotebooks\Scripts\python.exe C:\Users\hp\Documents\GitHub\Circuit_Notebooks\sample_obj_function.py

    str = strcat("C:\Users\hp\Documents\GitHub\CircuitNotebooks\Scripts\python.exe ", "C:\Users\hp\Documents\GitHub\Circuit_Notebooks\imfil_obj_function.py ");
    str = strcat(str, num2str(x(1), 12), " ", num2str(x(2), 12));
    disp(str);
    system(str);
    
    % Read in results
    outID = fopen("C:\Users\hp\Documents\GitHub\Circuit_Notebooks\imfil\output.txt",'r');
    [A,count] = fscanf(outID,"%f");
    
    disp(A(1));
    
    fout = A(1);
    ifail = A(2);
    icount = A(3);
    fclose(outID);
end
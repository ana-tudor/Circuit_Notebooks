import matlab.io.hdf4.*

format long
 
% This script is designed to bridge between imfil.m and sample_obj_function.py

% Set options if desired, have not seen them have an effect
options=imfil_optset;
% options=imfil_optset('termtol',10,options);
% options=imfil_optset('function_delta',-20,options);
% options=imfil_optset('verbose',1,options);
% options=imfil_optset('scaledepth',10,options);
% options=imfil_optset('scalestart',1,options);
% options=imfil_optset('fscale',-1.01,options);
% options=imfil_optset('random_stencil',2,options);
% options=imfil_optset('add_new_directions','tangent_directions',options);


% Read in initial values
ampID = fopen('/home/ana-tudor/Circuit_Notebooks/imfil/16_params_init.txt');
initamps = [];

% init = strrep(fgetl(ampID),';',' ');
init = fgetl(ampID);
disp(init)
while ischar(init)
    disp(init);
    initvec = reshape(str2num(char(strsplit(init))),16,[]);
    disp(initvec);
    initamps = [initamps initvec];
%     init = strrep(fgetl(ampID),';',' ');
    init = fgetl(ampID);
end

fclose(ampID);

disp(initamps)
outname='16_params_out.csv';

initamps2 = [];
results = [];
fvalues = [];
iterations = [];

%Try multiple initial values 
n=1;
while n < 11
    disp(size(initamps));
    initvec = initamps(1:end,n);
    disp(initvec);

    %Inner loop feeds results of previous minimization to check for
    %convergence
    inner=1;
    while inner < 4
    
        initamps2 = [initamps2; initvec.'];
        bounds1=[-1. 1.;-1. 1.;-1. 1.;-1. 1.]
        bounds2=[-2. 2.;-1. 1.;-2. 2.;-1. 1.]
        bounds=[bounds1;bounds2;bounds2;bounds2]
        [result, histout, comphist] = imfil(initvec, @func, 50, bounds, options);

        results = [results; result.'];
        fvalues = [fvalues; histout(end, 2)];
        iterations = [iterations; histout(end,1)];


        initvec = result;
        disp(initvec);
        disp(histout);
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

init_label = ["init1" "init2" "init3" "init4"];
init_label = [init_label "init5" "init6" "init7" "init8"]
init_label = [init_label "init9" "init10" "init11" "init12"]
init_label = [init_label "init13" "init14" "init15" "init16"]

results_label = ["fin1" "fin2" "fin3" "fin4"];
results_label = [results_label "fin5" "fin6" "fin7" "fin8"]
results_label = [results_label "fin9" "fin10" "fin11" "fin12"]
results_label = [results_label "fin13" "fin14" "fin15" "fin16"]

fid = fopen(outname, 'w') ;
fprintf(fid, '%s,', [init_label results_label "fvalues"]) ;
fprintf(fid, '%s\n', "iters") ;
fclose(fid) ;


dlmwrite(outname, [initamps2 results fvalues iterations], 'delimiter', ',', 'precision', 10, '-append'); 

% fclose(outID);

function [fout, ifail, icount] = func(x,h)
    % Run external python objective function here, through command line
    % integrated into matlab. Include x as cmd line arguments
    
    %!C:\Users\hp\Documents\GitHub\CircuitNotebooks\Scripts\python.exe C:\Users\hp\Documents\GitHub\Circuit_Notebooks\sample_obj_function.py

    str = strcat("/home/ana-tudor/quantum/bin/python3.6 ", "/home/ana-tudor/Circuit_Notebooks/imfil/16_param_obj_func.py ");
    str = strcat(str, num2str(x(1), 12), " ", num2str(x(2), 12), " ");
    str = strcat(str, num2str(x(3), 12), " ", num2str(x(4), 12), " ");
    str = strcat(str, num2str(x(5), 12), " ", num2str(x(6), 12), " ");
    str = strcat(str, num2str(x(7), 12), " ", num2str(x(8), 12), " ");
    str = strcat(str, num2str(x(9), 12), " ", num2str(x(10), 12), " ");
    str = strcat(str, num2str(x(11), 12), " ", num2str(x(12), 12), " ");
    str = strcat(str, num2str(x(13), 12), " ", num2str(x(14), 12), " ");
    str = strcat(str, num2str(x(15), 12), " ", num2str(x(16), 12));
    disp(str);
    system(str);
    
    % Read in results
    outID = fopen("/home/ana-tudor/Circuit_Notebooks/imfil/out_param_energies.txt",'r');
    [A,count] = fscanf(outID,"%f");
    
    disp(A(1));
    
    fout = A(1);
    ifail = A(2);
    icount = A(3);
    fclose(outID);
end
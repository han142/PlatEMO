%% Add path

addpath(genpath('MM_testfunctions/'));
addpath(genpath('IDMP_testfunctions/'));
addpath(genpath('Indicator_calculation/'));
addpath(genpath('MOEADBDN/'));
addpath(genpath('DNEA/'));
addpath(genpath('MMO_DE_CSCD/'));
addpath(genpath('MO_Ring_CSO_SCD/'));
addpath(genpath('MMODE_ICD/'));
addpath(genpath('NMOHSA/'));
addpath(genpath('NDNSGAII/'));
addpath(genpath('MMOPIO/'));
addpath(genpath('MMEA-SND/'));

clear all
clc

global alg_name fname

% 定义保存路径
base_save_dir = './Results'; % 基础保存目录
if ~exist(base_save_dir, 'dir')
    mkdir(base_save_dir); % 创建基础目录
end

N_function =  34 ;% number of test function
runtimes=1;  % odd number
local_function = [10,11,12,13,15,19]; % the indicator values of these functions are calculated using only global PS and PF
%% Initialize the parameters in MMO test functions
for i_func=8:8
    
    switch i_func
        case 1
            fname='MMF1';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[1 -1];     % the low bounds of the decision variables
            xu=[3 1];      % the up bounds of the decision variables
            repoint=[1.1,1.1]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 2
            fname='MMF2';
            n_obj=2;
            n_var=2;
            xl=[0 0];
            xu=[1 2];
            repoint=[1.1,1.1];
        case 3
            fname='MMF3';
            n_obj=2;
            n_var=2;
            xl=[0 0];
            xu=[1 1.5];
            repoint=[1.1,1.1];
        case 4
            fname='MMF4';
            n_obj=2;
            n_var=2;
            xl=[-1 0];
            xu=[1 2];
            repoint=[1.1,1.1];
        case 5
            fname='MMF5';
            n_obj=2;
            n_var=2;
            xl=[1 -1];
            xu=[3 3];
            repoint=[1.1,1.1];
        case 6
            fname='MMF6';
            n_obj=2;
            n_var=2;
            xl=[1 -1];
            xu=[3 2];
            repoint=[1.1,1.1];
        case 7
            fname='MMF7';
            n_obj=2;
            n_var=2;
            xl=[1 -1];
            xu=[3 1];
            repoint=[1.1,1.1];
        case 8
            fname='MMF8';
            n_obj=2;
            n_var=2;
            xl=[-pi 0];
            xu=[pi 9];
            repoint=[1.1,1.1];
        case 9
            fname='MMF9';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[0.1 0.1];     % the low bounds of the decision variables
            xu=[1.1 1.1];      % the up bounds of the decision variables
            repoint=[1.21,11]; % reference point used to calculate the Hypervolume
        case 10
            fname='MMF10';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[0.1 0.1];     % the low bounds of the decision variables
            xu=[1.1 1.1];      % the up bounds of the decision variables
            repoint=[1.21,13.2]; % reference point used to calculate the Hypervolume
        case 11
            fname='MMF11';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[0.1 0.1];     % the low bounds of the decision variables
            xu=[1.1 1.1];      % the up bounds of the decision variables
            repoint=[1.21,15.4];
        case 12
            fname='MMF12';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[0 0];     % the low bounds of the decision variables
            xu=[1 1];      % the up bounds of the decision variables
            repoint=[1.54,1.1];
        case 13
            %*need to be modified
            fname='MMF13';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[0.1 0.1 0.1];     % the low bounds of the decision variables
            xu=[1.1 1.1 1.1];      % the up bounds of the decision variables
            repoint=[1.54,15.4];
        case 14
            fname='MMF14';  % function name
            n_obj=3;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[0 0 0];     % the low bounds of the decision variables
            xu=[1 1 1];      % the up bounds of the decision variables
            repoint=[2.2,2.2,2.2];
        case 15
            fname='MMF15';  % function name
            n_obj=3;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[0 0 0];     % the low bounds of the decision variables
            xu=[1 1 1];      % the up bounds of the decision variables
            repoint=[2.5,2.5,2.5];
        case 16
            fname='MMF1_z';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[1 -1];     % the low bounds of the decision variables
            xu=[3 1];      % the up bounds of the decision variables
            repoint=[1.1,1.1];
        case 17
            fname='MMF1_e';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[1 -20];     % the low bounds of the decision variables
            xu=[3 20];      % the up bounds of the decision variables
            repoint=[1.1,1.1];
        case 18
            fname='MMF14_a';  % function name
            n_obj=3;
            n_var=3;
            xl=[0 0 0];
            xu=[1 1 1];
            repoint=[2.2,2.2,2.2];
        case 19
            fname='MMF15_a';  % function name
            n_obj=3;
            n_var=3;
            xl=[0 0 0];
            xu=[1 1 1];
            repoint=[2.5,2.5,2.5];
        case 20
            fname='SYM_PART_simple';
            n_obj=2;
            n_var=2;
            xl=[-20 -20];
            xu=[20 20];
            repoint=[4.4,4.4];
        case 21
            fname='SYM_PART_rotated';
            n_obj=2;
            n_var=2;
            xl=[-20 -20];
            xu=[20 20];
            repoint=[4.4,4.4];
        case 22
            fname='Omni_test';
            n_obj=2;
            n_var=3;
            xl=[0 0 0];
            xu=[6 6 6];
            repoint=[4.4,4.4];
        case 23
            fname='IDMPM2T1';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[-1,-1];     % the low bounds of the decision variables
            xu=[1,1];      % the up bounds of the decision variables
            repoint=[0.22,0.22]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 24
            fname='IDMPM2T2';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[-1,-1];     % the low bounds of the decision variables
            xu=[1,1];      % the up bounds of the decision variables
            repoint=[0.22,0.22]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 25
            fname='IDMPM2T3';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[-1,-1];     % the low bounds of the decision variables
            xu=[1,1];      % the up bounds of the decision variables
            repoint=[0.22,0.22]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 26
            fname='IDMPM2T4';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[-1,-1];     % the low bounds of the decision variables
            xu=[1,1];      % the up bounds of the decision variables
            repoint=[0.22,0.22]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 27
            fname='IDMPM3T1';  % function name
            n_obj=3;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[-1,-1,-1];     % the low bounds of the decision variables
            xu=[1,1,1];      % the up bounds of the decision variables
            repoint=[0.2,0.2,0.2]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 28
            fname='IDMPM3T2';  % function name
            n_obj=3;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[-1,-1,-1];     % the low bounds of the decision variables
            xu=[1,1,1];      % the up bounds of the decision variables
            repoint=[0.2,0.2,0.2]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 29
            fname='IDMPM3T3';  % function name
            n_obj=3;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[-1,-1,-1];     % the low bounds of the decision variables
            xu=[1,1,1];      % the up bounds of the decision variables
            repoint=[0.2,0.2,0.2]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 30
            fname='IDMPM3T4';  % function name
            n_obj=3;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[-1,-1,-1];     % the low bounds of the decision variables
            xu=[1,1,1];      % the up bounds of the decision variables
            repoint=[0.2,0.2,0.2]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 31
            fname='IDMPM4T1';  % function name
            n_obj=4;       % the dimensions of the decision space
            n_var=4;       % the dimensions of the objective space
            xl=[-1,-1,-1,-1];     % the low bounds of the decision variables
            xu=[1,1,1,1];      % the up bounds of the decision variables
            repoint=[0.22,0.22,0.22,0.22]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 32
            fname='IDMPM4T2';  % function name
            n_obj=4;       % the dimensions of the decision space
            n_var=4;       % the dimensions of the objective space
            xl=[-1,-1,-1,-1];     % the low bounds of the decision variables
            xu=[1,1,1,1];      % the up bounds of the decision variables
            repoint=[0.22,0.22,0.22,0.22]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 33
            fname='IDMPM4T3';  % function name
            n_obj=4;       % the dimensions of the decision space
            n_var=4;       % the dimensions of the objective space
            xl=[-1,-1,-1,-1];     % the low bounds of the decision variables
            xu=[1,1,1,1];      % the up bounds of the decision variables
            repoint=[0.22,0.22,0.22,0.22]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 34
            fname='IDMPM4T4';  % function name
            n_obj=4;       % the dimensions of the decision space
            n_var=4;       % the dimensions of the objective space
            xl=[-1,-1,-1,-1];     % the low bounds of the decision variables
            xu=[1,1,1,1];      % the up bounds of the decision variables
            repoint=[0.22,0.22,0.22,0.22]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
    
    end

    %% 选择优化算法并设置名称
    alg_name = 'MOEADBDN'; % 根据实际调用的算法动态设置

    % 创建算法子文件夹
    alg_save_dir = fullfile(base_save_dir, 'MOEADBDN_difference_size');
    if ~exist(alg_save_dir, 'dir')
        mkdir(alg_save_dir);
    end

    % 创建测试函数子文件夹
    func_save_dir = fullfile(alg_save_dir, fname, '100');
    if ~exist(func_save_dir, 'dir')
        mkdir(func_save_dir);
    end

    %% Load reference PS and PF data
    if i_func<=22
        if sum(i_func == local_function) >0
            load (strcat([fname,'_globalPSPF']))
            load (strcat([fname,'_localPSPF']));
            PS = global_PS;
            PF = global_PF;
        else
            load  (strcat([fname,'_Reference_PSPF_data']));
        end
    else
        load  (strcat(fname));
        PS = PSS;
    end
    %% Initialize the population size and the maximum evaluations
    if n_obj == 2
        popsize=100;
    elseif n_obj == 3
        popsize=406;
    elseif n_obj == 4
        popsize=560;
    end

    Max_fevs=120*1000;
    Max_Gen=fix(Max_fevs/popsize);
    Indicator=[];
    
    for j=1:runtimes
        %% Search the PSs using NCDE
        fprintf('Running optimization algorithm: %s, Running test function: %s, times = %d \n', alg_name,fname,j);
        % 根据 alg_name 调用相应的算法
        tic;
        switch alg_name
            case 'MOEADTN'
                [ps, pf] = MOEADTN(fname, xl, xu, n_obj, popsize, Max_Gen);
            case 'DNEA'
                [ps, pf] = DNEA(fname, xl, xu, n_obj, popsize, Max_Gen);
            case 'MMO_DE_CSCD'
                [ps, pf] = MMO_DE_CSCD(fname, xl, xu, n_obj, popsize, Max_Gen);
            case 'MO_Ring_CSO_SCD'
                [ps, pf] = MO_Ring_CSO_SCD(fname, xl, xu, n_obj, popsize, Max_Gen);
            case 'MMODE_ICD'
                [ps, pf] = MMODE_ICD(fname, xl, xu, n_obj, popsize, Max_Gen);
            case 'NMOHSA'
                [ps, pf] = NMOHSA(fname, xl, xu, n_obj, popsize, Max_Gen);
            case 'NDNSGAII'
                [ps, pf] = Decision_niched_NSGAII(fname, xl, xu, n_obj, popsize, Max_Gen);
            case 'MOEADBDN'
                [ps, pf] = MOEADBDN(fname, xl, xu, n_obj, popsize, Max_Gen);
            otherwise
                error('Unknown algorithm: %s', alg_name);
        end
        Time = toc; % 停止计时并获取时间
        fprintf('代码运行时间: %.3f 秒\n', Time);

        % Indicators
        hyp=Hypervolume_calculation(pf,repoint);
        IGDx=IGD_calculation(ps,PS);
        IGDf=IGD_calculation(pf,PF);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;%
        Indicator.(alg_name)(j, :) = [1./PSP, 1./hyp, IGDx, IGDf, Time];
        fprintf('Running test function: %s \n %d times rPSP=%f \n', fname, j, 1./PSP);
        PSdata.(alg_name){j} = ps;
        PFdata.(alg_name){j} = pf;
        clear ps pf hyp IGDx IGDf CR PSP Time
        
    end
    
    %% Choose one PS with the median indicators
    % Choose one PS for MMO_DE_CSCD
    choosen_In = Indicator.(alg_name)(:, 1); % Choose PS according to 1/PSP
    median_index = find(choosen_In == min(choosen_In));
    choose_ps.(alg_name){1} = PSdata.(alg_name){median_index};
    choose_pf.(alg_name){1} = PFdata.(alg_name){median_index};

    choosen_In = Indicator.(alg_name)(:, 2); % Choose PS according to 1/HV
    median_index = find(choosen_In == min(choosen_In));
    choose_ps.(alg_name){2} = PSdata.(alg_name){median_index};
    choose_pf.(alg_name){2} = PFdata.(alg_name){median_index};

    choosen_In = Indicator.(alg_name)(:, 3); % Choose PS according to IGDx
    median_index = find(choosen_In == min(choosen_In));
    choose_ps.(alg_name){3} = PSdata.(alg_name){median_index};
    choose_pf.(alg_name){3} = PFdata.(alg_name){median_index};

    choosen_In = Indicator.(alg_name)(:, 4); % Choose PS according to IGD
    median_index = find(choosen_In == min(choosen_In));
    choose_ps.(alg_name){4} = PSdata.(alg_name){median_index};
    choose_pf.(alg_name){4} = PFdata.(alg_name){median_index};

    choosen_In = Indicator.(alg_name)(:, 1)+Indicator.(alg_name)(:, 2)+Indicator.(alg_name)(:, 3)+Indicator.(alg_name)(:, 4); % Choose PS according to PSP
    median_index = find(choosen_In == min(choosen_In));
    choose_ps.(alg_name){5} = PSdata.(alg_name){median_index};
    choose_pf.(alg_name){5} = PFdata.(alg_name){median_index};

    clear choosen_In median_index
    
    %% Calculate mean and std of the indicators
    % NCDE
    Indicator.(alg_name)(runtimes+1, :) = min(Indicator.(alg_name)(1:runtimes, :)); %the minimum is the best
    Indicator.(alg_name)(runtimes+2, :) = max(Indicator.(alg_name)(1:runtimes, :)); %the maximum is the worst
    Indicator.(alg_name)(runtimes+3, :) = mean(Indicator.(alg_name)(1:runtimes, :));
    Indicator.(alg_name)(runtimes+4, :) = median(Indicator.(alg_name)(1:runtimes, :));
    Indicator.(alg_name)(runtimes+5, :) = std(Indicator.(alg_name)(1:runtimes, :));

    % Generate Table data in the report
    Table.(alg_name).rPSP(:, i_func) = Indicator.(alg_name)(runtimes+1:runtimes+5, 1);%Talbe II data
    Table.(alg_name).rHV(:, i_func) = Indicator.(alg_name)(runtimes+1:runtimes+5, 2);%Talbe III data
    Table.(alg_name).IGDX(:, i_func) = Indicator.(alg_name)(runtimes+1:runtimes+5, 3);%Talbe IV data
    Table.(alg_name).IGDF(:, i_func) = Indicator.(alg_name)(runtimes+1:runtimes+5, 4);%Talbe V data
    Table.(alg_name).IGDF(:, i_func) = Indicator.(alg_name)(runtimes+1:runtimes+5, 5);%Talbe V data

    %% save resultdata
    save(fullfile(func_save_dir, [fname, 'PSPF_indicator_global_data_', alg_name]), 'PSdata', 'PFdata', 'Indicator');
    save(fullfile(func_save_dir, [fname, 'ChoosenPSPFdata_', alg_name]), 'choose_ps', 'choose_pf');
    clear PSdata PFdata
    
end
save(fullfile(func_save_dir, [fname, 'Table_global_', alg_name]), 'Table');


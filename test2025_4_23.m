% 读取Excel文件
data = readtable('取对数加1.xls');

% 提取唯一年份
years = unique(data.year);
T = length(years); % 年份总数（2013到2023，共11年）

% 每个年份的DMU数量
n_dmu_per_year = 109; % 每个年份109个DMU
n_inputs = 2; % 两个投入：x1, x2
n_outputs = 1; % 一个产出：y

% 初始化Malmquist指数矩阵
% 109个DMU × (T-1)个年份对（2013-2014到2022-2023）
M = zeros(n_dmu_per_year, T-1);

% 技术参考集（所有年份的DMU）
X_ref = [data.x1, data.x2]; % 所有DMU的投入
Y_ref = data.y; % 所有DMU的产出

% 遍历连续年份对
for t = 1:(T-1)
    year_t = years(t);
    year_tp1 = years(t+1);
    
    % 提取年份 t 和 t+1 的数据
    data_t = data(data.year == year_t, :);
    data_tp1 = data(data.year == year_tp1, :);
    
    % 确保每个年份有109个DMU
    if height(data_t) ~= n_dmu_per_year || height(data_tp1) ~= n_dmu_per_year
        error('年份 %d 或 %d 的DMU数量不等于109', year_t, year_tp1);
    end
    
    % 按id排序以确保DMU对齐
    data_t = sortrows(data_t, 'id');
    data_tp1 = sortrows(data_tp1, 'id');
    
    % 遍历每个DMU
    for i = 1:n_dmu_per_year
        % 时间 t 和 t+1 的数据（DMU i）
        X_t = [data_t.x1(i), data_t.x2(i)]; % 年份 t 的投入
        Y_t = data_t.y(i); % 年份 t 的产出
        X_tp1 = [data_tp1.x1(i), data_tp1.x2(i)]; % 年份 t+1 的投入
        Y_tp1 = data_tp1.y(i); % 年份 t+1 的产出
        
        % 计算距离函数
        % D^t(x^t, y^t)：时间 t 在 t 技术下的效率
        Dt_xt_yt = dea_output_oriented(X_t, Y_t, X_ref, Y_ref);
        % D^{t+1}(x^{t+1}, y^{t+1})：时间 t+1 在 t+1 技术下的效率
        Dtp1_xtp1_ytp1 = dea_output_oriented(X_tp1, Y_tp1, X_ref, Y_ref);
        % D^t(x^{t+1}, y^{t+1})：时间 t+1 在 t 技术下的效率
        Dt_xtp1_ytp1 = dea_output_oriented(X_tp1, Y_tp1, X_ref, Y_ref);
        % D^{t+1}(x^t, y^t)：时间 t 在 t+1 技术下的效率
        Dtp1_xt_yt = dea_output_oriented(X_t, Y_t, X_ref, Y_ref);
        
        % 计算Malmquist指数（按照图片公式）
        % 第一部分：效率变化 (EC)
        EC = Dtp1_xtp1_ytp1 / Dt_xt_yt;
        % 第二部分：技术变化 (TC)
        TC = sqrt((Dt_xtp1_ytp1 / Dtp1_xtp1_ytp1) * (Dt_xt_yt / Dtp1_xt_yt));
        % Malmquist指数
        M(i, t) = EC * TC;
    end
end

% 创建表格（横坐标为年份对，纵坐标为id）
% 构造列名（2013-2014, 2014-2015, ...）
col_names = cell(1, T-1);
for t = 1:(T-1)
    col_names{t} = sprintf('%d-%d', years(t), years(t+1));
end

% 构造表格
results = array2table(M, 'VariableNames', col_names);
results.id = (1:n_dmu_per_year)';
results = movevars(results, 'id', 'Before', 1); % 将id列移到第一列

% 在命令行窗口显示表格（部分显示，避免输出太长）
disp('Malmquist指数表格（部分显示）：');
disp(results(1:5, :)); % 只显示前5行，避免输出过长

% 保存完整表格到Excel文件
writetable(results, 'Malmquist_Results_Table.xlsx');
fprintf('完整表格已保存到 Malmquist_Results_Table.xlsx\n');

% 产出导向DEA效率计算函数（VRS）
function theta = dea_output_oriented(Xo, Yo, X, Y)
    n = size(X, 1); % 决策单元数量
    m = size(X, 2); % 投入数量
    s = size(Y, 2); % 产出数量
    
    % 线性规划设置（最大化theta）
    f = [1; zeros(n, 1)]; % 目标：最大化theta
    Aeq = [zeros(m, 1), X'; -Yo, Y']; % 等式约束：投入和产出
    beq = [Xo'; 0]; % 等式约束右侧
    A = [0, ones(1, n)]; % VRS约束：lambda之和=1
    b = 1;
    lb = [0; zeros(n, 1)]; % 下界：theta >= 0, lambdas >= 0
    
    % 使用linprog求解线性规划
    options = optimoptions('linprog', 'Display', 'none');
    [x, fval] = linprog(-f, A, b, Aeq, beq, lb, [], options);
    
    theta = x(1); % 效率得分 (theta)
end

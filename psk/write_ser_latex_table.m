function latex_table_str = write_ser_latex_table(theo, sim, SNR, caption, label, filename)
%write_ser_latex_table 
%   Generates a latex table file for the SER values.

    latex_dir = 'latex/';
    fn = fullfile(latex_dir);
    if ~exist(fn, 'dir')
       mkdir(latex_dir);
    end

    % Caption and label.
    caption = strcat('\\caption{', caption ,'}');
    label = strcat('\\label{', label ,'}');
    
    % Percentage difference between theoretical and simulated
    diff = percentage_difference(theo, sim);
    
    % Build data rows for given SNR values.
    data_row = '';
    for index = SNR
        data_row = strcat(data_row, num2str(index), ' & ', num2str(theo(index+1)), ' & ', num2str(sim(index+1)), ' & ', num2str(diff(index+1)), '\\\\ \\hline');
    end
    
    % Rap data roes into table LaTex structure.
    latex_table_str = strcat('\\begin{table}[H]',...
        '\\centering',...
        caption,...
        label,...
        '\\begin{tabular}{|c|c|c|c|}',...
        '\\hline',...
        '\\multirow{2}{*}{\\textbf{SNR}} & \\multicolumn{3}{c|}{\\textbf{SER}} \\\\ \\cline{2-4}',... 
        '& \\textbf{Theoretical} & \\textbf{Simulated} & \\textbf{Difference} \\\\ \\hline',...
        data_row,...
        '\\end{tabular}',...
        '\\end{table}');
    
    % Write to .tex file
    f = fopen(strcat(latex_dir, filename), 'wt');
    fprintf(f, latex_table_str); % <- every three entries are written as a line
    fclose(f);
end


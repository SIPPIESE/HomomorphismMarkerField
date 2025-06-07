clear;
clc;


% This script file plots the demo coding images corresponding to window
% length n, demo-images of window length 5, 6, 7 are included in the demo
% coding directory ./demo_coding

% Here, select the available window length, recommended 4 - 7
n = 11;
select = [0, 0];

disp("Generating demo coding of window length [" + string(n) + "]");

% DATA READ
X = readmatrix('data/hmf/' + string(n) + '+' + string(n) + '_X.txt');
Y = readmatrix('data/hmf/' + string(n) + '+' + string(n) + '_Y.txt');
G = readmatrix('data/hmf/' + string(n) + '+' + string(n) + '_G.txt');
w = size(G, 1);
sample = n + 1;
min_x = select(1) - floor((sample - 1) / 2);
max_x = select(1) + (sample - floor((sample - 1) / 2) - 1);
min_y = select(2) - floor((sample - 1) / 2);
max_y = select(2) + (sample - floor((sample - 1) / 2) - 1);
flag_select = false;
if (min_x > 0 && max_x < w && min_y > 0 && max_y < w)
    flag_select = true;
end

% PLOT CONFIGURATION
font_size = 20;
figure('Position', [0, 0, 800, 800]);
hold on
set(gca, 'defaulttextinterpreter', 'latex');
set(gca, 'defaultAxesTickLabelInterpreter', 'latex');
set(gca, 'defaultLegendInterpreter', 'latex');

% DRAW RECTANGLE
explicit = '#fff';
implicit = '#000';
if (flag_select)
    explicit = '#ddd';
    implicit = '#ccc';
end
font_size = 20 * 12 / w;
for i = 1:w
    for j = 1:w
        if G(i, j) == 0
            rectangle('Pos', [i, - j - 1, 1, 1], 'FaceColor', implicit, 'EdgeColor', 'None');
            text(i + 0.5, -j - 0.53, '$\textbf{0}$', 'Color', explicit, 'FontSize', font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        else
            rectangle('Pos', [i, - j - 1, 1, 1], 'FaceColor', explicit, 'EdgeColor', 'None');
            text(i + 0.5, -j - 0.53, '$\textbf{1}$', 'Color', implicit, 'FontSize', font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        disp("Progress [" + string(roundn(((i - 1) * w + j) / (w * w) * 100, -1)) + "%]");
    end
end

% DRAW SELECT
if (flag_select)
    quiver(min_x - 1, -select(2)-0.5, sample + 2, 0, 1, 'Color', '#000', 'LineWidth', 1);
    quiver(select(1) + 0.5, -max_y - 2, 0, sample + 2, 1, 'Color', '#000', 'LineWidth', 1);
    for i = 1:w
        for j = 1:w
            if (i >= min_x && i <= max_x && j == select(2) || j >= min_y && j <= max_y && i == select(1))
                if G(i, j) == 0
                    rectangle('Pos', [i, - j - 1, 1, 1], 'FaceColor', '#000', 'EdgeColor', 'None');
                    text(i + 0.5, -j - 0.53, '$\textbf{0}$', 'Color', '#fff', 'FontSize', font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                else
                    rectangle('Pos', [i, - j - 1, 1, 1], 'FaceColor', '#fff', 'EdgeColor', 'None');
                    text(i + 0.5, -j - 0.53, '$\textbf{1}$', 'Color', '#000', 'FontSize', font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end
            end
        end
    end
end


% DRAW AXIS
line([1, w + 1], [-1, -1], 'Color', '#000', 'LineWidth', 1);
line([1, 1], [-1, - w - 1], 'Color', '#000', 'LineWidth', 1);
quiver(1, - w + 7, 0, -1, 8, 'Color', '#000', 'LineWidth', 1);
quiver(w - 7, -1, 1, 0, 8, 'Color', '#000', 'LineWidth', 1);

for j = 1:w
    if G(1, j) == 0
        text(0.5, -j - 0.53, '$\textbf{0}$', 'Color', '#000', 'FontSize', font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    else
        text(0.5, -j - 0.53, '$\textbf{1}$', 'Color', '#000', 'FontSize', font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
for i = 1:w
    if G(i, 1) == 0
        text(i + 0.5, - 0.53, '$\textbf{0}$', 'Color', '#000', 'FontSize', font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    else
        text(i + 0.5, - 0.53, '$\textbf{1}$', 'Color', '#000', 'FontSize', font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end


% PLOT CONFIGURATION
set(gca, 'xtick', [], 'ytick', [], 'xcolor', 'none', 'ycolor', 'none')
set(gca, 'position', [0, 0, 1, 1]);
axis equal;
axis([0, w + 1, - w - 1, 0]);

% EXPORT
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
filename = "plt/plt_" + string(n) + "+" + string(n) + ".pdf";
if (flag_select)
filename = "plt/plt_" + string(n) + "+" + string(n) + "_select.pdf";
end
print(filename, '-dpdf')
disp("Demo coding pattern saved!");
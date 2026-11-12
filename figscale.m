axs = findall(gcf, 'Type', 'axes');
for k = 1:numel(axs)
    if strcmp(axs(k).YLabel.String, 'Neuron ID')
        set(axs(k), 'CLim', [0 0.7]);
    end
end
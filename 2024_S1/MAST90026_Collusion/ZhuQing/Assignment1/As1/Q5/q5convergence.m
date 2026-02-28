N = 5:5:150;
figure; % Create a new figure
hold on; % Hold the plot to add multiple lines

for i = 1:numel(N) % Iterate over each value of N
    node = -1:1/N(i):1;
    s = length(node);
    [pt, qt, ft] = Assignment_1_Q5_BVP_Functions(node);
    [D, q, f] = toSelfAdjointForm(node, pt, qt, ft);
    elem = [linspace(1,s-1,s-1); linspace(2,s,s-1)]';
    u = BvpFE(node, elem, D, q, f, 2, 0);
    plot(node, u, 'DisplayName', ['N = ', num2str(N(i))]); % Set legend label
end

hold off; % Release the hold on the plot
legend('Location', 'best'); % Add legend
xlabel('x'); % Label for x-axis
ylabel('u(x)'); % Label for y-axis
title('Solution of Boundary Value Problem for Different Numbers of Nodes'); % Title for the plot

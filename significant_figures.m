function [rounded_numbers] = significant_figures(numbers)
%SIGNIFICANT_FIGURES calculates rounded numbers considering 2 significant
%figures of the minimum number
    precision=floor(log10(min(numbers)))-1;
    rounded_numbers=round(numbers/10^precision)*10^precision;
end


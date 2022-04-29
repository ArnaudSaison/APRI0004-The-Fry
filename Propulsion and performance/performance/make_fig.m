function make_fig(fig_name, Name, height, width, AR)

print_figures.H = height * 10;
print_figures.W = width * 10 * AR;
margin = 0;
name = [Name '.pdf'];

set(fig_name, 'PaperPosition', [margin/2 margin/2 print_figures.W print_figures.H])
set(fig_name, 'PaperSize', [print_figures.W+margin print_figures.H+margin])
%get(gca,'OuterPosition')
print(fig_name, name, '-dpdf', '-r300' );

end
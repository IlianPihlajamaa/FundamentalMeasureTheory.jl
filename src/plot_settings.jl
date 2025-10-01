mytheme = Theme(Axis=(xtickalign=1,
    xminortickalign=1,
    xticksmirrored=true,
    ytickalign=1,
    yminortickalign=1,
    yticksmirrored=true,
    xgridvisible=false,
    ygridvisible=false,))
set_theme!(merge(mytheme, theme_latexfonts()))
CairoMakie.theme(:palette).color[] = ColorSchemes.:ColorSchemes.:seaborn_muted6

colors = [ColorSchemes.:Dark2_8[(i-1)/7] for i in 1:8]
linestyles = [:solid, (:dash, :dense), :dash, :dashdot, (:dashdot, :dense), (:dot, :dense), :dot, (:dashdotdot, :dense)]

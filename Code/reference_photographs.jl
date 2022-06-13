using Plots

pgfplotsx()
plot_both = scatter([511], [510], xaxis=("x pixel coordinate",(1,1020)), yaxis=("y pixel coordinate",(1,1020)), label="Didymos", legend = :topright)
scatter!(plot_both, [937], [675], label="Dimorphos")

plot_didymos = scatter([780], [730], xaxis=("x pixel coordinate",(1,1020)), yaxis=("y pixel coordinate",(1,1020)), label="Didymos", legend = :topright)
scatter!(plot_didymos, [1937], [1675], label="Dimorphos")

plot_dimorphos = scatter([1780], [1730], xaxis=("x pixel coordinate",(1,1020)), yaxis=("y pixel coordinate",(1,1020)), label="Didymos", legend = :topright)
scatter!(plot_dimorphos, [372], [975], label="Dimorphos")

plot_none = scatter([1780], [1730], xaxis=("x pixel coordinate",(1,1020)), yaxis=("y pixel coordinate",(1,1020)), label="Didymos", legend = :topright)
scatter!(plot_none, [1372], [1975], label="Dimorphos")

savefig(plot_both, ".\\Results\\reference_both.pdf")
savefig(plot_didymos, ".\\Results\\reference_didymos.pdf")
savefig(plot_dimorphos, ".\\Results\\reference_dimorphos.pdf")
savefig(plot_none, ".\\Results\\reference_none.pdf")
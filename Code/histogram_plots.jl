using Plots, LaTeXStrings

pgfplotsx()

sma = [1180.339682,1180.278227,1180.372943,1180.182732,1180.368972,1180.483087,1180.434957,1180.385632,1180.381408,1180.376835,1180.573045,1180.353934,1180.352732,1180.443805,1180.447496,1180.299602,1180.327886,1180.299734,1180.355087,1180.356485]

ecc = [8.13E-06,1.04E-06,5.43E-06,8.27E-06,8.66E-06,8.64E-06,9.62E-06,1.95E-06,5.77E-07,2.72E-06,9.38E-06,3.28E-06,8.27E-06,8.71E-06,9.62E-06,4.12E-06,4.58E-06,8.28E-06,8.64E-06,2.70E-06]

lambda  = [147.305859,147.3281903,147.315742,147.3394979,147.35682,147.3136775,147.3167355,147.3407745,147.3149884,147.3142202,147.3313253,147.3386222,147.3394979,147.342159,147.3167355,147.3265379,147.3315979,147.3252727,147.3136775,147.3084908]

mu = [36.47655705,36.66749527,36.61636663,36.01666168,35.96385716,35.96385716,36.03227588,36.70674464,36.6960425,35.92487839,35.86587316,36.00496616,36.01666168,36.1027755,36.03227588,36.4683235,36.21157991,36.46385387,35.96385716,35.92487839]

hist_plot_sma = histogram(sma, bins=5, xlabel = "Semi-major axis, " * L"a,\ [m]", ylabel = "Frequency", legend=false)

hist_plot_ecc = histogram(ecc, bins=5, xlabel = "Eccentricity, " * L"e,\ [\ ]", ylabel = "Frequency", legend=false)

hist_plot_lambda = histogram(lambda, bins=5, xlabel = L"\textnormal{True Longitude},\ \lambda_{true},\ [\deg]", ylabel = "Frequency", legend=false)

hist_plot_mu = histogram(mu, bins=5, xlabel = L"\textnormal{Gravitational parameter},\ \mu,\ [m^{3}s^{-2}]", ylabel = "Frequency", legend=false)

savefig(hist_plot_sma, ".\\Results\\hist_plot_sma.pdf")
savefig(hist_plot_ecc, ".\\Results\\hist_plot_ecc.pdf")
savefig(hist_plot_lambda, ".\\Results\\hist_plot_lambda.pdf")
savefig(hist_plot_mu, ".\\Results\\hist_plot_mu.pdf")
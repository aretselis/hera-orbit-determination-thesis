using Plots, CSV, DataFrames, Dates, LaTeXStrings

function main()
    # Read file
    file_name = "C:\\Users\\retse\\repos\\hera-data\\Discovery_StatisticsPrintDownload.csv"
    df = DataFrame(CSV.File(file_name))
    end_date_index = 160
    ticks = minimum(df.Date[1:end_date_index]):Year(14):maximum(df.Date[1:end_date_index])
    pgfplotsx()
    image = plot(df.Date[1:end_date_index], df.NEA[1:end_date_index], xlims=(minimum(df.Date[1:end_date_index]), maximum(df.Date[1:end_date_index])), xticks = ticks, widen=true, formatter=:plain, legend = :topleft, label = "Near-Earth Asteroids")
    plot!(image, df.Date[1:end_date_index], df."NEA-km"[1:end_date_index], label = "Near-Earth Asteroids (over 1 km)")
    plot!(image, df.Date[1:end_date_index], df."NEA-140m"[1:end_date_index], label = "Near-Earth Asteroids (over 140 m)")
    plot!(image, df.Date[1:end_date_index], df.PHA[1:end_date_index], label= "Potential Hazardous Asteroids", xlabel = "Date", ylabel = "Number of objects")
    #display(image)
    savefig(".\\Results\\jpl_nea_data.pdf")
end

main()
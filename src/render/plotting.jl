function plot_position_overlap(ms, sort_perm, pics_olap_folder; width_factor=33)
    olap_ratio_ij = get_overlap_ratio(ms)
    olap_ratio_ij = olap_ratio_ij[sort_perm, sort_perm]
    x_ticks_1 = [i for i in 1:(ms.num_motifs)]; x_ticks_2 = ["D$i" for i in 1:ms.num_motifs]
    for i = 1:ms.num_motifs
        mask = (1:ms.num_motifs .!= i);
        x = collect(1:ms.num_motifs)[mask]
        y = @view olap_ratio_ij[i,:][mask]
        width = max(width_factor*ms.num_motifs, 400)
        fig = Figure(resolution = (width, 400), fonts = (; regular= "sans"));
        ax = Axis(fig[1, 1]; xticklabelrotation = pi/4, ylabel = "Jaccard index");
        ax.xticks = (x_ticks_1[mask], x_ticks_2[mask])
        ax.xticklabelsize = 18; ax.xlabelsize = 25; 
        ax.yticklabelsize = 18; ax.ylabelsize = 25;        
        ax.yticks = ([0.0, 0.25, 0.5, 0.75, 1.0], ["0.0", "0.25", "0.5", "0.75", "1.0"])
        ylims!(ax, high=1.0)
        barplot!(ax, x, y; strokewidth = 1, strokecolor = :black)
        save(joinpath(pics_olap_folder, "olap_d$i.png"), fig, px_per_unit = 0.6)
    end
end

function print_cmat_at_folder(cmats, folder; gt=false)
    char_ = gt ? "g" : "d"
    println("printing $(length(cmats)) count matrices at $folder")
    pfms = [cmats[i] ./ sum(cmats[i], dims=1) for i in 1:length(cmats)];
    logo_folder = folder
    for i in eachindex(pfms)
        io = open(logo_folder*"/d$i.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        g = Int.(floor.(pfms[i] .* 1000)); # make it a count matrix
        for j = 1:size(pfms[i],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(g[1,j])\t$(g[2,j])\t$(g[3,j])\t$(g[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i).transfac -n 40 --errorbars NO -F png --fineprint " " --resolution 72 -s large --fontsize 16 --color-scheme classic -o $(logo_folder)/$(char_)$(i).png`);
    end
end



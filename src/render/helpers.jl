make_grey(s::String) = grey_tag_front*s*grey_tag_back

get_folder_names(target_folder::String) = 
    target_folder*"/"*logo_olap_folder_name, 
    target_folder*"/"*logo_no_olap_folder_name,
    target_folder*"/"*pics_olap_folder_name,
    target_folder*"/"*pics_no_olap_folder_name

function make_folder_paths(folders::Vector{String})
    for folder in folders
        !isdir(folder) && mkpath(folder)
    end
end

function get_rounded_pval(pval::Real, low_pval)
    str = "$pval"; s = nothing;
    if !occursin("e-", str)
        s =  string(round(pval, sigdigits=3));        
    else
        q = split(str, "e-"); 
        q1len = length(q[1]);
        s = join([q[1][1:min(q1len, 4)], q[2]], "e-");
    end
    return !low_pval ? make_grey(s) : s;
end

function save_pfms_as_transfac(logo_folder::String, cmats, sort_perm::Vector{Int}, numbers)
    pfms = countmat2pfm.(cmats)

    @floop for (i,ind) in zip(numbers, sort_perm)
        this_pwm_size = "medium"
        io = open(logo_folder*"/d$i.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = Int.(floor.(pfms[ind] .* 100)); # make it a count matrix
        for j = 1:size(pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)

        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i).transfac -n 150 --number-fontsize 17 --errorbars NO -F png --fineprint " " --resolution 96 -s $this_pwm_size --fontsize 24 --small-fontsize 18 --color-scheme classic -o $(logo_folder)/d$(i).png`);

        # do it for the reverse complement as well
        io = open(logo_folder*"/d$(i)_c.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = Int.(floor.(reverse(pfms[ind]) .* 100));
        for j = 1:size(pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i)_c.transfac -n 150 --number-fontsize 17 --errorbars NO -F png --fineprint " " --resolution 96 -s $this_pwm_size --fontsize 24 --small-fontsize 18 --color-scheme classic -o $(logo_folder)/d$(i)_c.png`);
    end
end
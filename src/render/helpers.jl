make_grey(s::String) = grey_tag_front*s*grey_tag_back

get_folder_names(target_folder) = 
    target_folder*"/"*logo_olap_folder_name, 
    target_folder*"/"*logo_no_olap_folder_name,
    target_folder*"/"*pics_olap_folder_name,
    target_folder*"/"*pics_no_olap_folder_name

function make_folder_paths(folders)
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

function save_pfms_as_transfac(logo_folder::String, 
                               cmats, sort_perm::Vector{Int}, 
                               numbers;
                               data_name="unspecified")
    pfms = countmat2pfm.(cmats)
    counts_each_pfm = Int.(floor.([sum(c[:,1], dims=1)[1] for c in cmats]))

    @floop for (i,ind) in zip(numbers, sort_perm)
        this_pfm = pfms[ind];
        num_sites = counts_each_pfm[ind]
        this_pwm_size = "medium"
        pfm_len = size(this_pfm, 2)
        io = open(logo_folder*"/d$i.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = Int.(floor.(this_pfm .* num_sites)); # make it a count matrix
        for j = 1:pfm_len
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
        q = Int.(floor.(reverse(this_pfm) .* num_sites));
        for j = 1:pfm_len
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i)_c.transfac -n 150 --number-fontsize 17 --errorbars NO -F png --fineprint " " --resolution 96 -s $this_pwm_size --fontsize 24 --small-fontsize 18 --color-scheme classic -o $(logo_folder)/d$(i)_c.png`);

        # write as meme file as well
        io = open(logo_folder*"/d$i.meme", "w")
        print(io, "MEME version 4\n\n")
        print(io, "ALPHABET= ACGT\n\n")
        print(io, "strands: + -\n\n")
        print(io, "Background letter frequencies\n")
        print(io, "A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        print(io, "MOTIF $i $data_name \n")
        print(io, "letter-probability matrix: alength= 4 w= $pfm_len nsites= $num_sites E= 0\n")
        for i in axes(this_pfm,2)
            print(io, " ")
            for j = 1:4
                print(io, "$(this_pfm[j,i]) ")
            end
            print(io, "\n")
        end
        close(io)
    end
end
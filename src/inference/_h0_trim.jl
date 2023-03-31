function cmat2ic(cmat; bg=_bg_, ps = _pseudocounts_) # TODO change ps
    cmat_w_ps = cmat .+ ps
    freq_mat = cmat_w_ps ./ sum(cmat_w_ps, dims=1)
    reshape(sum(freq_mat .* log2.(freq_mat ./ bg), dims=1), size(cmat_w_ps,2))
end

function ic_span(ic_vec; ic_thresh=ic_trim_thresh, length_thresh=pfm_minimal_length)
    ic_vec_len = length(ic_vec)
    span_start = 1; 
    span_end = ic_vec_len
    for i = 1:ic_vec_len
        if ic_vec[i] < ic_thresh
            span_start += 1
        else
            break
        end
    end
    for i = ic_vec_len:-1:1
        if ic_vec[i] < ic_thresh
            span_end -= 1
        else
            break
        end
    end
    keep = span_end-span_start+1 ≥ length_thresh
    return span_start, span_end, keep
end

function trim_cmats(cmats, bg)
    new_cmats = Vector{eltype(cmats)}()
    ic = cmat2ic.(cmats; bg=bg)
    trim_info = ic_span.(ic)
    for (i, (span_start, span_end, keep)) in enumerate(trim_info)
        keep && push!(new_cmats, cmats[i][:,span_start:span_end])        
    end
    return new_cmats
end

function trim_H(data, H, bg)
    count_matrices = obtain_count_matrices(data, H)
    ic = cmat2ic.(count_matrices; bg=bg)
    trim_info = ic_span.(ic)
    new_dict = Dictionary{composition_key_type, Vector{value_type}}();
    for ((i, (span_start, span_end, keep)), k) in zip(enumerate(trim_info), keys(H))
        if keep 
            vs = Vector{value_type}(undef, length(H[k]))
            len = span_end - span_start + 1
            new_k = (f1=k.f1, f2=k.f2, f3=k.f3, d12=k.d12, d13=k.d13, len=len)
            for (ind,v) in enumerate(H[k])            
                vs[ind] = (seq_num=v.seq_num, pos=v.pos+span_start-1, comp=v.comp)
            end
            insert!(new_dict, new_k, vs)
        end
    end
    println("dict length : $(length(new_dict))")

    return new_dict
end

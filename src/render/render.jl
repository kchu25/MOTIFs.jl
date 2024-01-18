function render_main_summary_page_olap(labels, 
                                  pvalues, 
                                  logos,
                                  target_folder,
                                  target_folder_olap_pwms,
                                  target_folder_pics,
                                  valid_alphas,
                                  activate_counts,
                                  num_seqs
                                  )
    df = DataFrame(label=labels, eval=pvalues, logo=logos, counts=activate_counts);
    if length(valid_alphas) > 0
        out = Mustache.render(html_template_olap, 
                          target_folder=target_folder,
                          target_folder_pics=target_folder_pics,
                          valid_alphas="$valid_alphas",
                          num_alphas="$(length(valid_alphas)-1)",
                          min_alpha="$(valid_alphas[1])",
                          num_seq=num_seqs,
                          pic_folder=target_folder_pics,
                          logo_folder=target_folder_olap_pwms, DF=df);
    else
        out = Mustache.render(html_template_olap, 
                          target_folder=target_folder, num_seq=num_seqs,
                          pic_folder=target_folder_pics,
                          logo_folder=target_folder_olap_pwms, DF=df);
    end
    io = open(target_folder*"/summary.html", "w")
    print(io, out);
    close(io)
end

function render_main_summary_page_no_olap(labels, 
                                  pvalues, 
                                  logos,
                                  target_folder,
                                  target_folder_no_olap_pwms,
                                  valid_alphas,
                                  activate_counts,
                                  num_seqs
                                  )
    df = DataFrame(label=labels, eval=pvalues, logo=logos, counts=activate_counts);
    if length(valid_alphas) > 0
        out = Mustache.render(html_template_no_olap, 
                          target_folder=target_folder, 
                          valid_alphas="$valid_alphas",
                          num_alphas="$(length(valid_alphas)-1)",
                          min_alpha="$(valid_alphas[1])",
                          num_seq=num_seqs,
                          logo_folder=logo_folder_name, DF=df);
    else
        out = Mustache.render(html_template_no_olap, 
                          target_folder=target_folder, num_seq=num_seqs,
                          logo_folder=target_folder_no_olap_pwms, DF=df);
    end
    io = open(target_folder*"/summary_no_olap.html", "w")
    print(io, out);
    close(io)
end

function render_result!(target_folder, ms, data, bg; alpha_fisher = 1e-5)

    logo_olap_folder, logo_no_olap_folder, pics_olap_folder, pics_no_olap_folder = get_folder_names(target_folder);
    make_folder_paths([target_folder, logo_olap_folder, logo_no_olap_folder, pics_olap_folder, pics_no_olap_folder]);

    order_for_display = obtain_groupings_for_display1(ms);

    ############## overlap version render ##############
    @info "Scanning the foreground..."
    scan_w_gpu!(ms, data);
    @info "Scanning the shuffled background..."
    scan_w_gpu!(ms, data; bg=true);
    @time filter_positions_scores_usecomp!(ms, data, bg);
    active_counts, _ = 
        get_uniq_counts(ms)

    @info "Calculating p-values..."
    pvalues, sort_perm_1, uniq_active_counts_test = 
        get_pvec_and_related_info2(ms, data, alpha_fisher, order_for_display)
    labels  = ["D$j" for j = 1:ms.num_motifs];
    logo_names   = ["d$(j)" for j = 1:ms.num_motifs];
    # @info "plot the overlap..."
    @info "save the PWMs..."
    save_pfms_as_transfac(logo_olap_folder, ms.cmats, sort_perm_1, collect(1:ms.num_motifs));
    activate_counts_total = (Int.(active_counts .+ uniq_active_counts_test))[sort_perm_1]
    valid_alphas = Int[]; # TODO: verify the purpose of this and remove it if necessary
 
    render_main_summary_page_olap(labels, 
                                pvalues, 
                                logo_names,
                                target_folder,
                                logo_olap_folder_name,
                                pics_olap_folder_name,
                                valid_alphas,
                                activate_counts_total,
                                data.N+data.N_test
                                );
end

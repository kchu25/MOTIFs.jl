function discover_motifs(datapath, save_path;
                         num_epochs=nothing)
    @info "load data"
    data = FASTA_DNA{float_type}(datapath)
    this_bg = get_data_bg(data)
    @info "training..."
    cdl, hp, len, projs = train_ucdl(data; num_epochs=num_epochs)
    @info "extract motifs..."
    ms = run_thru(data, cdl, hp, len, projs, this_bg)
    render_result!(save_path, ms, data, this_bg)
end
function discover_motifs(datapath, save_path;
                         num_epochs=nothing)
    data = FASTA_DNA{float_type}(datapath)
    this_bg = get_data_bg(data)
    cdl, hp, len, projs = train_ucdl(data; num_epochs=num_epochs)
    ms = run_thru(data, cdl, hp, len, projs, this_bg)
    render_result!(save_path, ms, data, this_bg)
end
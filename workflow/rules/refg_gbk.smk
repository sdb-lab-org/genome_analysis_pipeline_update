def get_ref_paths(wildcards):
    ref_dir = os.path.join(WORKING_DIR, "reference")
    standard_gbk = os.path.join(ref_dir, "refg.gbk")
    standard_fasta = os.path.join(ref_dir, "refg.fasta")

    if os.path.exists(standard_gbk) and os.path.exists(standard_fasta):
        print(f"✔ Using manually provided reference: {standard_gbk}, {standard_fasta}")
        return standard_gbk, standard_fasta

    gbk_candidates = [f for f in os.listdir(ref_dir) if f.endswith(".gbk")]
    fasta_candidates = [f for f in os.listdir(ref_dir) if f.endswith(".fasta") or f.endswith(".fa") or f.endswith(".fna")]

    if gbk_candidates and fasta_candidates:
        gbk = os.path.join(ref_dir, gbk_candidates[0])
        fasta = os.path.join(ref_dir, fasta_candidates[0])
        print(f"✔ Found existing files: {gbk}, {fasta}")
        return gbk, fasta

    ref_list = config.get("refg", [])
    ref_sample = ref_list[0] if ref_list else sorted(ALL_GENOMES)[0]
    gbk = os.path.join(OUTPUT_DIR, ref_sample, "annotation", "prokka", "prokka.gbk")
    fasta = os.path.join(OUTPUT_DIR, ref_sample, "annotation", "prokka", "prokka.fna")
    print(f"✔ Using Prokka output from {ref_sample}")
    return gbk, fasta


rule refg_gbk:
    input:
        get_ref_paths
    output:
        gbk = WORKING_DIR + "reference/refg.gbk",
        fasta = WORKING_DIR + "reference/refg.fasta"
    message:
        "Copying reference .gbk and .fasta to standard locations"
    run:
        import shutil
        import os

        gbk_in, fasta_in = input
        if os.path.abspath(gbk_in) != os.path.abspath(output.gbk):
            print(f"✔ Copying {gbk_in} → {output.gbk}")
            shutil.copy(gbk_in, output.gbk)
        else:
            print("✔ refg.gbk already in place — no copy needed.")

        if os.path.abspath(fasta_in) != os.path.abspath(output.fasta):
            print(f"✔ Copying {fasta_in} → {output.fasta}")
            shutil.copy(fasta_in, output.fasta)
        else:
            print("✔ refg.fasta already in place — no copy needed.")    

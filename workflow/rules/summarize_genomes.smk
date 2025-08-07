rule summarize_sample_stats:
    input:
        quast_reports      = expand(OUTPUT_DIR + "01_qc/quast/{sample}/report.txt", sample=SAMPLES),
        annotation_reports = expand(OUTPUT_DIR + "08_temp/temp_txt/{sample}.txt", sample=SAMPLES)
    output:
        summary = OUTPUT_DIR + "01_qc/summary_qc.tsv"
    message:
        "Creating summary with GC, contig count, total length and CDS for each sample (Prokka or Bakta)"
    run:
        import os
        import re
        import pandas as pd

        summary_data = []

        # Map sample → quast report
        quast_map = {}
        for path in input.quast_reports:
            match = re.search(r'/quast/([^/]+)/report\.txt', path)
            if match and os.path.exists(path) and os.path.getsize(path) > 0:
                quast_map[match.group(1)] = path
            else:
                print(f"[WARN] QUAST fehlt oder leer: {path}")

        # Map sample → annotation report (could be Prokka or Bakta)
        anno_map = {}
        for path in input.annotation_reports:
            match = re.search(r'/temp_txt/([^/]+)\.txt$', path)
            if match and os.path.exists(path) and os.path.getsize(path) > 0:
                anno_map[match.group(1)] = path
            else:
                print(f"[WARN] Annotation fehlt oder leer: {path}")

        valid_samples = sorted(set(quast_map) & set(anno_map))

        if not valid_samples:
            raise ValueError("No valid QUAST or annotation reports found")

        for sample in valid_samples:
            row = {"sample_name": sample}

            # ---- Parse QUAST ----
            try:
                with open(quast_map[sample]) as f:
                    for line in f:
                        if "GC (%)" in line:
                            row["GC (%)"] = line.strip().split()[-1]
                        elif "# contigs" in line:
                            row["# contigs"] = line.strip().split()[-1]
                        elif "Total length" in line:
                            row["Total length"] = line.strip().split()[-1]
            except Exception as e:
                print(f"[ERROR] QUAST parsing error for {sample}: {e}")
                continue

            # ---- Parse Annotation (Prokka or Bakta) ----
            try:
                with open(anno_map[sample]) as f:
                    content = f.read()

                if "CDS:" in content and "contigs:" in content:
                    # Prokka format
                    for line in content.splitlines():
                        if line.startswith("CDS:"):
                            row["CDS"] = line.strip().split()[-1]

                elif "CDSs:" in content and "Sequence(s):" in content:
                    # Bakta format
                    for line in content.splitlines():
                        if line.startswith("CDSs:"):
                            row["CDS"] = line.strip().split()[-1]
                        elif line.startswith("GC:") and "GC (%)" not in row:
                            row["GC (%)"] = line.strip().split()[-1]
                else:
                    print(f"[WARN] Unknown annotation format for {sample}")

            except Exception as e:
                print(f"[ERROR] Annotation parsing error for {sample}: {e}")
                continue

            summary_data.append(row)

        df = pd.DataFrame(summary_data)
        df.to_csv(output.summary, sep="\t", index=False)

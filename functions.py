import re


def parse_NPDC_DIAMOND_blastp_tabular(blastp_tabular_file):
    blastp_results = []
    with open(blastp_tabular_file, "r", encoding="utf-8") as f:
        for l in f:
            if l.startswith("#"):
                match = re.match(r"^# [0-9]+ hits found(.*)", l)
                if not match:
                    continue
                else:
                    l = match.group(1).strip()
            parts = l.strip().split("\t")
            if len(parts) < 9:
                continue
            blastp_results.append(
                {
                    "query_id": parts[0],
                    "subject_id": int(parts[1]),  # CDS id of NPDC
                    "percent_identity": float(parts[2]),
                    "e_value": float(parts[3]),
                    "bit_score": float(parts[4]),
                    "query_start": int(parts[5]),
                    "query_end": int(parts[6]),
                    "subject_start": int(parts[7]),
                    "subject_end": int(parts[8]),
                }
            )
    return blastp_results


def get_target_info(target_nr, target_name, target_dir):
    print(f"  {target_nr}: {target_name}")
    blastp_meta = target_dir / f"metadata_{target_name}.tsv"
    blastp_tabular = target_dir / "blast_tabular_result.txt"
    assert (
        blastp_meta.is_file()
    ), f"Expected metadata file not found: {blastp_meta}"
    assert (
        blastp_tabular.is_file()
    ), f"Expected blast tabular file not found: {blastp_tabular}"
    return parse_NPDC_DIAMOND_blastp_tabular(blastp_tabular), blastp_meta


def cal_maxgap(locations):
    target_nrs = set(loc[4] for loc in locations)
    contigs = set(loc[0] for loc in locations)
    max_gaps = {}
    for contig in contigs:
        locations_in_contig = sorted(
            [loc for loc in locations if loc[0] == contig], key=lambda x: x[1]
        )  # Filter by contig
        targets_in_contig = set(loc[4] for loc in locations_in_contig)
        if len(targets_in_contig) < len(target_nrs):
            max_gaps[contig] = float(
                "inf"
            )  # Infinite gap if not all targets are present
            continue
        max_gap = 0
        for i in range(1, len(locations_in_contig)):
            gap = (
                locations_in_contig[i][1] - locations_in_contig[i - 1][2]
            )  # start of current - end of previous
            if gap > max_gap:
                max_gap = gap
        max_gaps[contig] = max_gap
    return max_gaps

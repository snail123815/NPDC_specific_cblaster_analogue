import re
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path

from functions import cal_maxgap, get_target_info

HELP_DESCRIPTION = """\
Process NPDC-specific cblaster results.

To prepare:
1. Run BLASTP for each target against the NPDC database. It allows
    5 targets per run maximum. Please use simple naming for targets as
    fasta headers (e.g., >TacA, >TacB, etc.)
    You cannot change NPDC DIAMOND-BLASTP parameters:
        e-value - 1e-10
        identity - 40%
        query coverage - 80%
2. Download the 'BLAST tabular result (.txt)' for every target. 
    You need to click download at least twice, first time clicking, the 
    server will prepare the data. Data can be only downloaded from a 
    second attempt. Expect many 'cooldown' timers between downloads.
3. Extract downloaded zip files, organize folders like the following,
    make sure the dir name per protein result is:
    NPDC + "_" + number + "_" + ProteinName(exactly from fasta header)

    SCRIPT_FOLDER/
    ├── NPDC_01_ProA/
    │   ├── blast_tabular_result.txt
    │   └── metadata_ProA.txt
    ├── NPDC_02_ProB/
    │   ├── blast_tabular_result.txt
    │   └── metadata_ProB.txt
    ├── NPDC_03_ProC/
    │   ├── blast_tabular_result.txt
    │   └── metadata_ProC.txt
    ├── Others.../
    ├── npdc_specific_cblaster_main.py
    └── functions.py

    The above naming convension matches the default --target_folder_pattern
    parsing rule: "^NPDC_([0-9]+)_(.+)$"

    Number should not be important, but it makes life easier to follow your
    input BGC ;-)

4. Run this script with appropriate arguments, e.g.:

    % python npdc_specific_cblaster_main.py \\
        --required_proteins TacA TacB TacC # the ProteinName from fasta.
        
    Output will be shown as text in your terminal.

5. Optionally download the genomes manually and do a real cblaster run to 
    generate a plot.
"""

parser = ArgumentParser(
    description=HELP_DESCRIPTION,
    formatter_class=RawDescriptionHelpFormatter,
)
parser.add_argument(
    "--target_folder_pattern",
    type=str,
    default=r"^NPDC_([0-9]+)_(.+)$",
    help=(
        "Regex pattern to identify target folders and extract target "
        "number and name. Default: '^NPDC_([0-9]+)_(.+)$' "
        "(matches folders like 'NPDC_1_TacA' and extracts target number '1' "
        "and target name 'TacA')."
    ),
)
parser.add_argument(
    "--required_proteins",
    nargs="*",
    default=[],
    help=(
        "List of required target protein names (e.g., TacA TacB TacC)."
        " Must match the target_name in the folder names (e.g., NPDC_1_TacA)."
    ),
)
parser.add_argument(
    "--gap",
    type=int,
    default=20000,
    help=(
        "Maximum allowed gap (in bp) between hits of different targets "
        "on the same contig. Default: 20000)"
    ),
)

args = parser.parse_args()
target_folder_regex = re.compile(args.target_folder_pattern)
required_proteins = args.required_proteins
filters = {
    "gap": args.gap,  # Maximum allowed gap (in bp) between hits of different targets on the same contig
    "min_unique_NI": 3,  # not implemented
    "min_identity_NI": 40.0,  # not implemented, must >= 40%
    "max_evalue_NI": 1e-10,  # not implemented, must >= 1e-10
    "min_qcov_NI": 80.0,  # not implemented, must >= 80%
}

op_targets = []
req_targets = []
for f in Path(".").iterdir():
    match = target_folder_regex.match(f.name)
    if f.is_dir() and match:
        target_nr, target_name = match.groups()
        target_nr = int(target_nr)
        if target_name in required_proteins:
            req_targets.append((target_nr, target_name, f))
        else:
            op_targets.append((target_nr, target_name, f))
op_targets.sort(key=lambda x: x[0])
req_targets.sort(key=lambda x: x[0])


print(f"Found {len(req_targets)} required targets:")
blastp_tabular_data = {}
for target_nr, target_name, target_dir in req_targets:
    blastp_tabular_data[target_name] = get_target_info(
        target_nr, target_name, target_dir
    )
print(f"Found {len(op_targets)} targets:")
for target_nr, target_name, target_dir in op_targets:
    blastp_tabular_data[target_name] = get_target_info(
        target_nr, target_name, target_dir
    )

genomes = {}
species_names = {}
for target_nr, target_name, target_dir in req_targets:
    tabular_data, blastp_meta = blastp_tabular_data[target_name]
    assert (
        tabular_data[0]["query_id"] == target_name
    ), f"Expected query ID {target_name} in blastp results"

    with open(blastp_meta, "r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue  # Skip header
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            bgc_nr = int(parts[0])
            genome_nr = int(parts[1])
            cds_nr = int(parts[2])
            contig_nr = int(parts[7])
            species_name = parts[4]
            strain_nr = int(parts[6])
            bgc_region = parts[11].split(":")[1]
            for blastp_hit in tabular_data:
                if blastp_hit["subject_id"] == cds_nr:
                    pidentity = blastp_hit["percent_identity"]
                    break
            else:
                raise ValueError(
                    f"No blastp hit found for CDS ID {cds_nr} in target {target_name}"
                )
            location = (
                contig_nr,
                int(parts[8]),
                int(parts[9]),
                int(parts[10]),
                target_nr,  # target number, not target name, for easier gap calculation
            )  # contig_nr, start, end, strand, query id
            ######## data strusture!!! #########
            target_info: dict = {
                "target_nr": target_nr,
                "target_name": target_name,
                "target_hits": {
                    "cds_nrs": [
                        cds_nr
                    ],  # Subject ID from blastp (CDS ID of NPDC)
                    "percent_identity": [pidentity],
                    "bgc_nrs": [bgc_nr],
                    "bgc_regions": [bgc_region],
                    "locations": [
                        location,
                    ],
                },
            }
            if genome_nr not in genomes:
                genomes[genome_nr] = [target_info]
                species_names[genome_nr] = [species_name, strain_nr]
            ###################################
            else:
                for existing_target in genomes[genome_nr]:
                    if existing_target["target_nr"] == target_nr:
                        for k, v in (
                            ("cds_nrs", cds_nr),
                            ("percent_identity", pidentity),
                            ("bgc_nrs", bgc_nr),
                            ("bgc_regions", bgc_region),
                            ("locations", location),
                        ):
                            existing_target["target_hits"][k].append(v)
                        break
                else:
                    genomes[genome_nr].append(target_info)
print(f"Total genomes with at least one required target hit: {len(genomes)}")

target_genomes = {}
print(
    f"""Filtering genomes
    1. with hits for all {len(req_targets)} required targets
    2. with less than {filters['gap']} bp gap between hits of different targets
    --- not implemented yet ---
    3. with at least {filters['min_unique_NI']} unique hits for each target (NI for not implemented)
    4. with at least {filters['min_identity_NI']}% identity for each hit (NI for not implemented)
    5. with at most {filters['max_evalue_NI']} e-value for each hit (NI for not implemented)
    6. with at least {filters['min_qcov_NI']}% query coverage for each hit (NI for not implemented)
"""
)


for k, v in genomes.items():
    # Check if all required targets are present
    if len(v) < len(req_targets):
        continue

    # Check gap between hits of different targets
    locations = []
    for target_info in v:
        for loc in target_info["target_hits"]["locations"]:
            locations.append(loc)
    max_gaps = cal_maxgap(locations)
    if not min(max_gaps.values()) < filters["gap"]:
        continue
    gap_on_contigs = {}
    for contig, gap in max_gaps.items():
        if gap >= filters["gap"]:
            continue
        else:
            gap_on_contigs[contig] = gap

    # Output genomes that pass all filters
    target_genomes[k] = (v, gap_on_contigs)

print()
print(f"Genomes with hits for all required targets: {len(target_genomes)}")
for genome_nr, (targets, gap_on_contigs) in target_genomes.items():
    print()
    print(f"Genome ID: {genome_nr}")
    print(f"  Species: {species_names[genome_nr][0]}")
    print(f"  Strain ID: {species_names[genome_nr][1]}")
    bgc_regions_total = set()
    for target_info in targets:
        for bgc_region in target_info["target_hits"]["bgc_regions"]:
            bgc_regions_total.add(bgc_region)
    print(f"  BGC Regions of hits: {', '.join(map(str, bgc_regions_total))}")
    print(
        f"  Contigs with gap < {filters['gap']} bp: {', '.join(map(str, gap_on_contigs.keys()))}"
    )


for target_nr, target_name, target_dir in op_targets:
    # with open(blastp_meta, "r") as f:
    pass

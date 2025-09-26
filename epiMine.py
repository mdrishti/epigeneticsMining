import requests
import pandas as pd
import argparse
import json
from typing import List, Dict, Any
from pathlib import Path
from io import StringIO
import re
from xml.etree import ElementTree as ET


ENA_API = "https://www.ebi.ac.uk/ena/portal/api/search"
EUROPEPMC_API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
FULLTEXT_API = "https://www.ebi.ac.uk/europepmc/webservices/rest/{id}/fullTextXML"


def query_ena(tax_id: int, platform: str) -> pd.DataFrame:
    # query ENA for read_experiment results by tax_tree and platform
    query = f"tax_tree({tax_id}) AND instrument_platform=\"{platform}\""
    data = {
        "result": "read_experiment",
        "query": query,
        #"fields": "experiment_accession,run_accession,study_accession,sample_accession,submission_accession,instrument_model",
        "fields": "age,aligned,altitude,assembly_quality,assembly_software,bam_aspera,bam_bytes,bam_ftp,bam_galaxy,bam_md5,base_count,binning_software,bio_material,bisulfite_protocol,broad_scale_environmental_context,broker_name,cage_protocol,cell_line,cell_type,center_name,checklist,chip_ab_provider,chip_protocol,chip_target,collected_by,collection_date,collection_date_end,collection_date_start,completeness_score,contamination_score,control_experiment,country,cultivar,culture_collection,datahub,depth,description,dev_stage,disease,dnase_protocol,ecotype,elevation,environment_biome,environment_feature,environment_material,environmental_medium,environmental_sample,experiment_accession,experiment_alias,experiment_target,experiment_title,experimental_factor,experimental_protocol,extraction_protocol,faang_library_selection,fastq_aspera,fastq_bytes,fastq_ftp,fastq_galaxy,fastq_md5,file_location,first_created,first_public,germline,hi_c_protocol,host,host_body_site,host_genotype,host_gravidity,host_growth_conditions,host_phenotype,host_scientific_name,host_sex,host_status,host_tax_id,identified_by,instrument_model,instrument_platform,investigation_type,isolate,isolation_source,last_updated,lat,library_construction_protocol,library_gen_protocol,library_layout,library_max_fragment_size,library_min_fragment_size,library_name,library_pcr_isolation_protocol,library_prep_date,library_prep_date_format,library_prep_latitude,library_prep_location,library_prep_longitude,library_selection,library_source,library_strategy,local_environmental_context,location,location_end,location_start,lon,marine_region,mating_type,ncbi_reporting_standard,nominal_length,nominal_sdev,pcr_isolation_protocol,ph,project_name,protocol_label,read_count,read_strand,restriction_enzyme,restriction_enzyme_target_sequence,restriction_site,rna_integrity_num,rna_prep_3_protocol,rna_prep_5_protocol,rna_purity_230_ratio,rna_purity_280_ratio,rt_prep_protocol,run_accession,run_alias,run_date,salinity,sample_accession,sample_alias,sample_capture_status,sample_collection,sample_description,sample_material,sample_prep_interval,sample_prep_interval_units,sample_storage,sample_storage_processing,sample_title,sampling_campaign,sampling_platform,sampling_site,scientific_name,secondary_project,secondary_sample_accession,secondary_study_accession,sequencing_date,sequencing_date_format,sequencing_location,sequencing_longitude,sequencing_method,sequencing_primer_catalog,sequencing_primer_lot,sequencing_primer_provider,serotype,serovar,sex,specimen_voucher,sra_aspera,sra_bytes,sra_ftp,sra_galaxy,sra_md5,status,strain,study_accession,study_alias,study_title,sub_species,sub_strain,submission_accession,submission_tool,submitted_aspera,submitted_bytes,submitted_format,submitted_ftp,submitted_galaxy,submitted_host_sex,submitted_md5,submitted_read_type,surveillance_target,tag,target_gene,tax_id,tax_lineage,taxonomic_classification,taxonomic_identity_marker,temperature,tissue_lib,tissue_type,transposase_protocol,variety",
        "format": "tsv",
    }
    resp = requests.post(ENA_API, data=data, timeout=1000)
    resp.raise_for_status()
    if not resp.text.strip():
        return pd.DataFrame()
    return pd.read_csv(StringIO(resp.text), sep="\t")


def query_europepmc(prj_id: str, tax_id: int, output_dir: Path) -> Dict[str, Any]:
    # query EuropePMC for a given PRJEB accession
    query = f"ACCESSION_TYPE:bioproject AND ACCESSION_ID:{prj_id}"
    params = {"query": query, "format": "json"}
    resp = requests.get(EUROPEPMC_API, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()
    has_hits = False
    results = []
    for record in data.get("resultList", {}).get("result", []):
        entry = record.get("pmcid")
        #    "pubYear": record.get("pubYear"),
        #}
        results.append(entry)
        
        # fetch full-text and mine flow cell info
        xml = fetch_fulltext(entry)
        if xml:
            has_hits = True
            hits = mine_flowcell(xml, entry, tax_id, output_dir)
            print(f"{prj_id}:{entry} Flow cell terms found: {hits}")
        else:
            print(f"{prj_id}:{entry} No OA full text available")

    return data, has_hits


def download_ena_runs_from_df(
    ena_df: pd.DataFrame,
    projects: List[str],
    tax_id: int,
    #genus_taxid: int,
    output_dir: Path
) -> None:
    # Download ENA run files for projects with flowcell info,
    # reusing the ENA dataframe already fetched in query_ena.
    for prj in projects:
        subset = ena_df[ena_df["study_accession"] == prj]
        if subset.empty:
            continue

        prj_dir = output_dir / f"tax_{tax_id}" / prj
        prj_dir.mkdir(parents=True, exist_ok=True)

        for _, row in subset.iterrows():
            run_acc = row["run_accession"]
            print(run_acc)
            ftp_field = row.get("fastq_ftp") or row.get("submitted_ftp")
            if pd.isna(ftp_field):
                continue
            urls = str(ftp_field).split(";")

            for url in urls:
                file_name = url.split("/")[-1]
                dest = prj_dir / file_name
                if dest.exists():
                    continue
                #ftp_url = f"https://ftp.sra.ebi.ac.uk/vol1/{url}"
                ftp_url = f"https://{url}"
                print(f"[INFO] Downloading {ftp_url} → {dest}")
                try:
                    r = requests.get(ftp_url, stream=True, timeout=600)
                    r.raise_for_status()
                    with open(dest, "wb") as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            f.write(chunk)
                except Exception as e:
                    print(f"[ERROR] Failed {ftp_url}: {e}")



def run_pipeline(tax_ids: List[int], platforms: List[str], output_dir: Path) -> None:
    # run ENA + EuropePMC pipeline for multiple tax IDs and platforms
    output_dir.mkdir(parents=True, exist_ok=True)
    europepmc_results = {}
    flowcell_positive_projects = []

    for tax_id in tax_ids:
        for platform in platforms:
            print(f"[INFO] Processing tax_id={tax_id}, platform={platform}")
            ena_df = query_ena(tax_id, platform)
            if ena_df.empty:
                print(f"[WARN] No ENA hits for tax_id={tax_id}, platform={platform}")
                continue

            ena_out = output_dir / f"ena_results_tax{tax_id}_{platform}.tsv"
            ena_df.to_csv(ena_out, sep="\t", index=False)
            print(f"[INFO] Saved ENA results → {ena_out}")
            cols = [c for c in ["study_accession", "tax_id", "tax_lineage"] if c in ena_df.columns]
            prj_ids = ena_df.loc[:, cols].dropna().drop_duplicates()
            #prj_ids = ena_df[["study_accession", "tax_id", "tax_lineage"]].dropna().drop_duplicates()
            #prj_ids = ena_df["study_accession"].dropna().unique()
            #for prj_id in prj_ids:
            for _, row in prj_ids.iterrows():
                prj_id = row["study_accession"]
                sample_tax_id = int(row["tax_id"])
                lineage = row["tax_lineage"]

                # derive genus taxid
                #genus_taxid = extract_genus_taxid(lineage, sample_tax_id)

                if prj_id not in europepmc_results:
                    europepmc_results[prj_id], has_hits = query_europepmc(prj_id, tax_id, output_dir)
                    if has_hits:
                        flowcell_positive_projects.append((prj_id, tax_id))
            for prj_id, tax_id in flowcell_positive_projects:
                download_ena_runs_from_df(ena_df, [prj_id], tax_id, output_dir)

    with open(output_dir / "europepmc_results.json", "w", encoding="utf-8") as f:
        json.dump(europepmc_results, f, indent=2)
    print(f"[INFO] Saved EuropePMC results → {output_dir/'europepmc_results.json'}")


def fetch_fulltext(pmid_or_pmcid: str) -> str:
    # fetch full text XML from EuropePMC if available
    url = FULLTEXT_API.format(id=pmid_or_pmcid)
    resp = requests.get(url, timeout=60)
    if resp.status_code == 200 and resp.text.strip():
        return resp.text
    return ""

def mine_flowcell(xml_text: str, entry: str, tax_id: int, output_dir: Path) -> List[str]:
    # extract possible flow cell numbers or IDs from full text XML
    #print(xml_text)
    tax_dir = output_dir / f"tax_{tax_id}"
    tax_dir.mkdir(parents=True, exist_ok=True)
    xml_path = tax_dir / f"{entry}.xml"
    #with open(f"{entry}.xml", "w") as f:
    with open(xml_path, "w", encoding="utf-8") as f:
        f.write(xml_text)

    tree = ET.fromstring(xml_text)
    text_chunks = [" ".join(elem.itertext()) for elem in tree.iter() if elem.text]
    text = " ".join(text_chunks)
    pattern = re.compile(r"FLO-\w+|R\d{1,2}\.\d{1,2}|MIN\d+|flow\s*cell", re.IGNORECASE) # regex for flowcells
    hits = list(set(pattern.findall(text)))
    tsv_path = tax_dir / "flowcell_hits.tsv"
    header = not tsv_path.exists()
    with open(tsv_path, "a", encoding="utf-8") as tsv_file:
        if header:
            tsv_file.write("pmc_id\tgenus_taxid\thits\n")
        tsv_file.write(f"{entry}\t{tax_id}\t{','.join(hits) if hits else 'NA'}\n")
    return hits





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate ENA + EuropePMC queries")
    parser.add_argument("--tax-ids", type=str, required=True, help="Comma-separated tax IDs or path to a file containing one ID per line")
    parser.add_argument("--platform", type=int, required = True, help="Choose platform set: 1=OXFORD_NANOPORE, 2=PACBIO_SMRT")
    parser.add_argument("--output-dir", type=str, default="results", help="Directory to save outputs")
    args = parser.parse_args()

    # parse tax_ids argument
    if Path(args.tax_ids).is_file():
        with open(args.tax_ids) as f:
            tax_ids = [int(line.strip()) for line in f if line.strip()]
    else:
        tax_ids = [int(x) for x in args.tax_ids.split(",") if x.strip()]

    # platforms to process
    platform_map = {1: ["OXFORD_NANOPORE"], 2: ["PACBIO_SMRT"]}
    platforms = platform_map[args.platform]

    run_pipeline(tax_ids, platforms, Path(args.output_dir))


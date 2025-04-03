# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "duckdb",
#     "requests",
# ]
# ///

"""Create test data from GWAS Catalog gene based submissions

Manually curated pubmed IDs:

34662886
36088354
36596879
37592023
39180217
39362880
40021682
40073867
36450978
37949852
39385933* (window-based)

and a bunch of GCST accessions (gcsts.txt)

Test data contains roughly 5% of the original files
"""

import requests
import duckdb
import re
import pathlib
import time

PUBMED_IDS = [
    34662886,
    36088354,
    36596879,
    37592023,
    39180217,
    39362880,
    40021682,
    40073867,
    36450978,
    37949852,
]


def fetch_study_links(pubmed_id):
    """Get study links recursively from a pubmed ID

    By default, only grab the first 100 GCSTs
    """
    BASE_URL = (
        "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId"
    )
    HEADERS = {"Accept": "application/json"}

    page = 0
    size = 50
    study_links = []

    while True:
        print(f"Fetching {pubmed_id=} {page=}...")

        # Send the request
        response = requests.get(
            BASE_URL,
            params={"pubmedId": pubmed_id, "page": page, "size": size},
            headers=HEADERS,
        )
        response.raise_for_status()
        data = response.json()

        # Extract matching studies
        studies = data.get("_embedded", {}).get("studies", [])
        for study in studies:
            disease_trait = study.get("diseaseTrait", {}).get("trait", "").lower()
            # magic strings that flag gene based sumstats in the trait field
            if (
                "gene-based" in disease_trait
                or "gene burden" in disease_trait
            ):
                study_links.append(study["_links"]["self"]["href"])

        # Check if there are more pages
        total_pages = data.get("page", {}).get("totalPages", 1)

        if page >= total_pages - 1:
            break

        page += 1

    return study_links


def get_directory_range(n):
    """Create an FTP directory range

    90083565 -> 90083001, 90084000
    90083000 -> 90082001, 90083000
    """
    if n % 1000 == 0:
        n -= 1
    lower = (n // 1000) * 1000 + 1
    upper = (n // 1000 + 1) * 1000
    return lower, upper


def make_ftp_url(gcst_url):
    """Make an FTP url from a GCST accession:"""
    base_url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

    # Extract GCST ID (e.g., "GCST90082112")
    match = re.search(r"GCST\d+", gcst_url)
    if not match:
        raise ValueError(f"GCST URL does not contain GCST {gcst_url}")

    gcst_id = match.group()

    # Determine the range (e.g., "GCST90082001-GCST90083000")
    gcst_num = int(gcst_id[4:])
    lower_bound, upper_bound = get_directory_range(gcst_num)

    range_folder = f"GCST{lower_bound:08d}-GCST{upper_bound:08d}"

    # Construct full URL
    full_url = f"{base_url}/{range_folder}/{gcst_id}"

    return full_url


def get_tsv_url(url):
    """Get a TSV file path (possibly gzip compressed) from a HTML directory listing

    >>> url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081596"
    >>> get_tsv_url(url)
    http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081596/GCST90081596_buildGRCh38.tsv.gz
    """
    response = requests.get(url)
    response.raise_for_status()

    # Use regex to find all TSV links
    try:
        tsv_file = re.findall(r'href="([^"]+\.(?:tsv\.gz|tsv))"', response.text)[0]
    except IndexError:
        raise ValueError(f"Couldn't find TSV file at {url}")

    time.sleep(0.2)  # naive rate limit for EBI load balancer
    return url + "/" + tsv_file


def sample_csv(path, outf, seed=42, overwrite=False):
    """Reproducibly sample 5% of rows from a text file and write to outf"""
    if pathlib.Path(outf).exists() and not overwrite:
        print(f"File {outf} already exists, skipping")
    else:
        print(f"Writing {outf}")
        with duckdb.connect(":memory:") as conn:
            sql = f"""
            COPY (SELECT * FROM read_csv(\"{path}\") 
            USING SAMPLE 5% (system, 42))        
            TO \"{outf}\"
            DELIMITER '\t';        
            """
            conn.execute(sql)


def get_sumstat_from_pubmed(pubmed_id):
    links = fetch_study_links(pubmed_id)
    test_data_path = pathlib.Path(__file__).parent.parent.resolve() / "tests" / "data"
    (test_data_path / str(pubmed_id)).mkdir(exist_ok=True, parents=True)

    for link in links:
        url = make_ftp_url(link)
        tsv_url = get_tsv_url(url)
        outf = str(test_data_path / str(pubmed_id) / f"{link.split("/")[-1]}.tsv.gz")
        sample_csv(path=tsv_url, outf=outf)

    return len(links)


def read_gcsts():
    with open("gcsts.txt") as f:
        next(f)  # skip the header comment
        return f.readlines()


def get_sumstat_from_gcst(gcst):
    test_data_path = (
        pathlib.Path(__file__).parent.parent.resolve() / "tests" / "data" / "gcsts"
    )
    test_data_path.mkdir(exist_ok=True, parents=True)
    url = make_ftp_url(gcst)
    tsv_url = get_tsv_url(url)
    outf = str(test_data_path / f"{gcst.split('/')[-1].strip()}.tsv.gz")
    sample_csv(path=tsv_url, outf=outf)


def main() -> None:
    for pubmed_id in PUBMED_IDS:
        n_processed = get_sumstat_from_pubmed(pubmed_id)
        if n_processed == 0:
            raise ValueError(f"No sumstats downloaded for {pubmed_id=}")

    gcsts = read_gcsts()
    for gcst in gcsts:
        get_sumstat_from_gcst(gcst)

if __name__ == "__main__":
    main()

# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "duckdb",
#     "requests",
#     "tenacity"
# ]
# ///

"""Create test data from GWAS Catalog gene based submissions

Using a manually curated list, query the Catalog API with pubmed IDs and a bunch of GCST accessions (gcsts.txt).

- Test data are gzip compressed TSV files
- These files contain roughly 5% of the original sumstat files
- Written to the tests/data directory

Run this script with uv to automatically install the dependencies:

$ uv run create_test_data.py
"""

import requests
import duckdb
import re
import pathlib
from tenacity import retry, wait_exponential

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


def get_gcst_from_pubmed(pubmed_id):
    """Get a list of gene-based GCSTs from a pubmed id"""
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

        studies = data.get("_embedded", {}).get("studies", [])
        for study in studies:
            disease_trait = study.get("diseaseTrait", {}).get("trait", "").lower()
            # magic strings that flag gene based sumstats in the trait field
            if "gene-based" in disease_trait or "gene burden" in disease_trait:
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


@retry(wait=wait_exponential(multiplier=1, min=2, max=10))  # (EBI load balancer)
def make_ftp_url(gcst_url):
    """Make an FTP url that links to a sumstats directory from a GCST accession"""
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


@retry(wait=wait_exponential(multiplier=1, min=2, max=10))  # (EBI load balancer)
def get_tsv_url(url):
    """Get a TSV URL (possibly gzip compressed) from an FTP directory listing"""
    response = requests.get(url)
    response.raise_for_status()

    # Use regex to find the first TSV link (there should only be one)
    try:
        tsv_file = re.findall(r'href="([^"]+\.(?:tsv\.gz|tsv))"', response.text)[0]
    except IndexError:
        raise ValueError(f"Couldn't find TSV file at {url}")

    return url + "/" + tsv_file


def sample_csv(path, outf, overwrite=False):
    """Reproducibly sample 5% of rows from a text file and write to a TSV"""
    if pathlib.Path(outf).exists() and not overwrite:
        print(f"File {outf} already exists, skipping")
    else:
        print(f"Writing {outf}")
        with duckdb.connect(":memory:") as conn:
            # duckdb fetches the remote TSV and decompresses it magically
            # seed = 42 for reproducibility
            sql = f"""
            COPY (SELECT * FROM read_csv(\"{path}\") 
            USING SAMPLE 5% (system, 42))        
            TO \"{outf}\"
            DELIMITER '\t';        
            """
            conn.execute(sql)


def get_sumstats_from_pubmed_id(pubmed_id):
    """Query the Catalog API with a pubmed ID and download sampled sumstats files"""
    links = get_gcst_from_pubmed(pubmed_id)
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
    """Query the Catalog API with a GCST accession and download sampled sumstats files"""
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
        n_processed = get_sumstats_from_pubmed_id(pubmed_id)
        if n_processed == 0:
            # no GCSTs were found from a pubmed ID. that's bad!
            raise ValueError(f"No sumstats downloaded for {pubmed_id=}")

    gcsts = read_gcsts()
    for gcst in gcsts:
        get_sumstat_from_gcst(gcst)


if __name__ == "__main__":
    main()

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
- These files contain 1000 rows from the original sumstat files
- Sampled data is written to the tests/data directory

Run this script with uv to automatically install the dependencies:

$ uv run create_test_data.py

or if using nox in the main package directory:

$ nox -s create_test_data
"""

import concurrent.futures
import gzip
import json
import os
import tempfile

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


def get_sumstats_from_pubmed_id(pubmed_id):
    """Query the Catalog API with a pubmed ID and return the sumstat TSV URL"""
    links = get_gcst_from_pubmed(pubmed_id)
    test_data_path = pathlib.Path(__file__).parent.parent.resolve() / "tests" / "data"
    (test_data_path / str(pubmed_id)).mkdir(exist_ok=True, parents=True)

    urls = [make_ftp_url(link) for link in links]
    return [get_tsv_url(url) for url in urls]


def get_sumstat_from_gcst(gcst):
    """Query the Catalog API with a GCST accession return the sumstat TSV URL"""
    url = make_ftp_url(gcst)
    return get_tsv_url(url)


def is_sumstat_ok(outf):
    """Sometimes sampled sumstat files get written with just a header :("""
    with gzip.open(outf, "rt") as f:
        n_lines = sum(1 for _ in f)
        if n_lines < 2:
            failed = True  # just a header or an empty file
        else:
            failed = False
    return failed


def sample_csv(path, outf, overwrite=False):
    """Reproducibly sample rows from a text file and write to a TSV"""
    outf_path = pathlib.Path(outf)

    if outf_path.exists() and not overwrite:
        print(f"File {outf} already exists, skipping")
    else:
        sql = f"""
        COPY (SELECT * FROM read_csv(\"{path}\", sample_size = 100_000) 
        USING SAMPLE reservoir(1000 ROWS) REPEATABLE (42))
        TO \"{outf}\"
        DELIMITER '\t';        
        """
        # only one memory database can exist a time
        # need to get a temporary file name for the duckdb database that doesn't exist
        fd, db_path = tempfile.mkstemp(prefix=outf_path.stem)
        os.close(fd)
        pathlib.Path(db_path).unlink()

        with duckdb.connect(db_path) as conn:
            # duckdb fetches the remote TSV and decompresses it magically
            # seed = 42 for reproducibility
            conn.execute(sql)

        pathlib.Path(db_path).unlink()

    failed = is_sumstat_ok(outf)

    if failed:
        print(f"Deleting bad sumstat file {path} and retrying")
        pathlib.Path(outf).unlink()
        raise ValueError(f"Empty output: {outf}")


def download_sumstats(tsv_urls):
    """Download and sample sumstats files in parallel"""
    print("Downloading sumstats...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for url_dict in tsv_urls:
            # make the directory for the data
            test_data_path = (
                pathlib.Path(__file__).parent.parent.resolve() / "tests" / "data"
            )
            gwas_id = url_dict["id"]
            (test_data_path / gwas_id).mkdir(exist_ok=True, parents=True)

            for url in url_dict["urls"]:
                basename = url.split("/")[-1].strip().split(".")[0]
                outf = str(test_data_path / gwas_id / f"{basename}.tsv.gz")
                futures.append(
                    executor.submit(sample_csv, path=url, outf=outf, overwrite=False)
                )

        for future in concurrent.futures.as_completed(futures):
            _ = future.result()  # correctly raise exceptions


def query_gwascatalog_api(gcst_path, tsv_path):
    """Get GCSTs from the GWAS Catalog API"""
    tsv_urls = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        with open(gcst_path) as f:
            next(f)  # skip the header comment
            gcsts = f.readlines()
        print("Querying GCSTs...")
        futures = [executor.submit(get_sumstat_from_gcst, gcst) for gcst in gcsts]
        for future in concurrent.futures.as_completed(futures):
            url = {"id": "gcsts", "urls": [future.result()]}
            tsv_urls.append(url)

    # get gcsts for pubmed IDs
    pubmed_links = [
        {"id": str(pubmed_id), "urls": get_sumstats_from_pubmed_id(pubmed_id)}
        for pubmed_id in PUBMED_IDS
    ]
    tsv_urls.extend(pubmed_links)

    with open(tsv_path, "w") as f:
        json.dump(tsv_urls, f, ensure_ascii=False, indent=4)


def main():
    script_dir = pathlib.Path(__file__).parent.resolve()
    gcsts = script_dir / pathlib.Path("gcsts.txt")
    tsv_urls = script_dir / pathlib.Path("urls.json")

    # get sumstat file URLs from the GWAS Catalog API
    if not tsv_urls.exists():
        query_gwascatalog_api(gcst_path=gcsts, tsv_path=tsv_urls)
    else:
        with open(tsv_urls) as f:
            tsv_urls = json.load(f)

    # randomly sample the URLs with duckdb and write sumstats to the tests data directory
    download_sumstats(tsv_urls)


if __name__ == "__main__":
    main()

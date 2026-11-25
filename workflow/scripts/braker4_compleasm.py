"""Prepare a run-private Compleasm cache and verify BRAKER4 rescue outputs."""

import configparser
import hashlib
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


def complete_summary(path):
    text = Path(path).read_text()
    for key in ("S", "D", "F", "M"):
        if not re.search(rf"(?m)^\s*{key}:\s*\d", text):
            raise ValueError(f"Incomplete Compleasm summary: {path}")
    return text


def private_cache(run_dir):
    return Path(run_dir) / "compleasm_downloads"


def container_command(args, run_dir):
    ini = configparser.ConfigParser()
    ini.read(args.template)
    uri = ini["containers"]["braker3_image"]
    image = Path(args.container_cache) / (hashlib.md5(uri.encode()).hexdigest() + ".simg")
    if not image.is_file():
        raise FileNotFoundError(f"Run the BRAKER4 container download first: {image}")
    return ["apptainer", "exec", "--cleanenv", "--no-home", "--bind", f"{run_dir},{args.container_cache}",
            str(image), "python3"]


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def verify_lineage(source, candidate):
    profiles = {p.name for p in (source / "hmms").glob("*.hmm")}
    if not profiles or profiles != {p.name for p in (candidate / "hmms").iterdir()}:
        raise ValueError(f"Compleasm and BUSCO HMM profile inventories differ: {candidate}")
    paths = ["scores_cutoff", "refseq_db.faa.gz", *(f"hmms/{name}" for name in sorted(profiles))]
    # The pinned Compleasm downloader only extracts the reference proteins,
    # HMMs and score cutoffs. BRAKER4 separately reads lengths_cutoff.
    if (candidate / "lengths_cutoff").exists():
        paths.append("lengths_cutoff")
    for name in paths:
        if digest(source / name) != digest(candidate / name):
            raise ValueError(f"Compleasm lineage differs from the selected BUSCO dataset: {name}")


def download_cache(args, run_dir, command):
    cache = private_cache(run_dir)
    source = Path(args.busco_downloads) / "lineages" / args.lineage
    if cache.exists():
        if list(cache.glob("*.tmp")):
            raise ValueError(f"Incomplete Compleasm download requires review: {cache}")
        verify_lineage(source, cache / args.lineage)
        if not (cache / args.lineage / "lengths_cutoff").is_file():
            raise FileNotFoundError(f"BRAKER4 lengths_cutoff is missing from {cache}")
        return cache
    with tempfile.TemporaryDirectory(dir=run_dir, prefix=".compleasm-cache-") as tmp:
        staged = Path(tmp) / "library"
        with (Path(run_dir) / "compleasm_download.log").open("w") as log:
            subprocess.run([*command, "/opt/compleasm_kit/compleasm.py", "download",
                            args.lineage, "-L", str(staged)], check=True, stdout=log,
                           stderr=subprocess.STDOUT, timeout=3600)
        lineage = staged / args.lineage
        verify_lineage(source, lineage)
        if not (lineage / "lengths_cutoff").exists():
            shutil.copyfile(source / "lengths_cutoff", lineage / "lengths_cutoff")
        staged.rename(cache)
    return cache


def prepare_cache(args, run_dir):
    """Use the pinned container's unmodified downloader before admitting a run."""
    command = container_command(args, run_dir)
    cache = download_cache(args, run_dir, command)
    # Exercise the same protein-mode startup used by best_by_compleasm.
    reference = cache / args.lineage / "refseq_db.faa.gz"
    probe = Path(run_dir) / "compleasm_preflight.faa"
    extract = (
        "import gzip,sys;"
        "f=gzip.open(sys.argv[1],'rt');o=open(sys.argv[2],'w');"
        "header=next(f);assert header.startswith('>');o.write(header);"
        "\nfor line in f:\n"
        " if line.startswith('>'): break\n"
        " o.write(line)\n"
        "o.close()"
    )
    subprocess.run([*command, "-c", extract, str(reference), str(probe)], check=True)
    with tempfile.TemporaryDirectory(dir=run_dir, prefix=".compleasm-probe-") as tmp:
        output = Path(tmp) / "results"
        with (Path(run_dir) / "compleasm_preflight.log").open("w") as log:
            subprocess.run([*command, "/opt/compleasm_kit/compleasm.py", "protein",
                            "-p", str(probe), "-l", args.lineage, "-L", str(cache),
                            "-o", str(output), "-t", "1"],
                           check=True, stdout=log, stderr=subprocess.STDOUT, timeout=600)
        complete_summary(output / "summary.txt")


def validate_rescue(run_dir, sample):
    root = Path(run_dir) / "output" / sample
    log = (root / "best_by_compleasm.log").read_text()
    if re.search(r"Return code of subprocess was [1-9]|Traceback|script may have failed", log):
        raise ValueError("BRAKER4 best_by_compleasm failed; exported predictions are not accepted")
    # Upstream suppresses subprocess failures. A result for every compared set
    # is required even when the original TSEBRA prediction remains the best.
    evidence = {}
    for name in ("genemark", "augustus", "braker"):
        folder = root / "best_by_compleasm_tmp" / "single" / name
        summaries = list(folder.rglob("summary.txt"))
        if len(summaries) != 1:
            raise ValueError(f"Missing Compleasm rescue evaluation for {name}: {folder}")
        evidence[name] = complete_summary(summaries[0])
    complete_summary(root / "compleasm_proteins" / "summary.txt")
    return evidence

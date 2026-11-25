"""Prepare a run-private Compleasm cache and verify BRAKER4 rescue outputs."""

import configparser
import hashlib
import re
import subprocess
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


def prepare_cache(args, run_dir):
    """Use the pinned container's unmodified downloader before admitting a run."""
    cache = private_cache(run_dir)
    cache.mkdir(exist_ok=True)
    if list(cache.glob("*.tmp")):
        raise ValueError(f"Incomplete Compleasm download requires review: {cache}")
    command = container_command(args, run_dir)
    with (Path(run_dir) / "compleasm_download.log").open("w") as log:
        subprocess.run([*command, "/opt/compleasm_kit/compleasm.py", "download",
                        args.lineage, "-L", str(cache)],
                       check=True, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
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
    output = Path(run_dir) / "compleasm_preflight"
    with (Path(run_dir) / "compleasm_preflight.log").open("w") as log:
        subprocess.run([*command, "/opt/compleasm_kit/compleasm.py", "protein",
                        "-p", str(probe), "-l", args.lineage, "-L", str(cache),
                        "-o", str(output), "-t", "1"],
                       check=True, stdout=log, stderr=subprocess.STDOUT, timeout=600)
    summaries = list(output.rglob("summary.txt"))
    if len(summaries) != 1:
        raise ValueError("Compleasm protein preflight did not produce one summary")
    complete_summary(summaries[0])


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

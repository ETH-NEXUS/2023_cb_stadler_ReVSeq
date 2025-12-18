#!/usr/bin/env python

import argparse
import sys
from io import StringIO
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate an EMBL flat file for one sample from a multi-FASTA "
            "consensus file, by fetching each segment's EMBL record from ENA "
            "and replacing ID + SQ."
        )
    )
    p.add_argument(
        "--consensus",
        required=True,
        help="Consensus FASTA for one sample (one record per segment; "
             "each header must contain the accession, e.g. KU509707.1).",
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output EMBL file containing all segments for that sample.",
    )
    p.add_argument(
        "--sample-id",
        required=True,
        help="Sample ID, e.g. 4L3pXZ (used as prefix in each ID line).",
    )
    return p.parse_args()

def fetch_from_ncbi_genbank(accession: str) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gb",
        "retmode": "text",
    }
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(params)
    sys.stderr.write(f"[INFO] Fetching GenBank from NCBI: {url}\n")
    with urlopen(url) as resp:
        text = resp.read().decode("utf-8")
    if not text.strip():
        raise RuntimeError(f"Empty GenBank response for accession {accession}")
    record = SeqIO.read(StringIO(text), "genbank")

    #print (text)  # DEBUG
    return record.format("embl")

def fetch_embl_text(accession: str) -> str:
    if accession.startswith("NC_"):
        return fetch_from_ncbi_genbank(accession)
    url = f"https://www.ebi.ac.uk/ena/browser/api/embl/{accession}"
    sys.stderr.write(f"[INFO] Fetching EMBL from ENA: {url}\n")
    try:
        with urlopen(url) as resp:
            text = resp.read().decode("utf-8")
    except HTTPError as e:
        raise RuntimeError(f"HTTP error fetching EMBL for {accession}: {e}")
    except URLError as e:
        raise RuntimeError(f"URL error fetching EMBL for {accession}: {e}")

    if not text.strip():
        raise RuntimeError(f"Empty EMBL response for accession {accession}")

    return text


def format_sq_block(seq: str):
    """
    Given a raw sequence string, return a list of EMBL SQ lines:
      - header 'SQ   Sequence ...'
      - wrapped sequence lines with 60 bp/line, groups of 10, final position right-aligned.
    """
    seq = "".join(seq.split())  # remove whitespace/newlines
    seq = seq.lower()           # ENA usually uses lowercase in SQ

    length = len(seq)
    a = seq.count("a")
    c = seq.count("c")
    g = seq.count("g")
    t = seq.count("t")
    other = length - (a + c + g + t)

    lines = []
    header = f"SQ   Sequence {length} BP; {a} A; {c} C; {g} G; {t} T; {other} other;"
    lines.append(header)

    # 60 bases per line, grouped into chunks of 10
    for i in range(0, length, 60):
        chunk = seq[i : i + 60]
        seqlen=len(chunk)
        if seqlen < 60: # pad with spaces if needed
            chunk = chunk.ljust(60)
        groups = [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
        body = " ".join(groups)
        pos = i + seqlen
        # 5 leading spaces, then sequence groups, then position right-aligned in width 9
        line = f"     {body} {pos:>9}"
        lines.append(line)

    return '\n'.join(lines)


def patch_single_embl(embl_text: str, sample_id: str, accession: str, consensus_seq: str) -> str:
    """
    Take the EMBL text for a single accession, and:
      - replace ID name with '<sample_id>_<accession>'
      - replace SQ section with one built from consensus_seq
      - ensure organism looks like influenza (string check).
    Return the patched EMBL text for that one record.
    """
    if "influenza a" not in embl_text.lower():
        raise RuntimeError(
            f"EMBL record for {accession} does not look like influenza."
        )

    lines = embl_text.rstrip("\n").split("\n")
    keep_tags = ['ID', 'AC', 'DE', 'OS', 'OC', 'FH', 'FT','XX','SQ','//']
    ft_keep={'source': ['organism', 'strain',  'serotype', 'segment']}
    current_ft=''
    keep_lines=[]
    sq_added=False
    for l in lines:
        tag = l[:2]
        if tag not in keep_tags:
            continue
        if tag =='ID':
            old_id=l[5:].strip()
            keep_lines.append(f'ID   {sample_id}_{old_id}')
        elif tag == 'FT':
            key = l[5:21].strip()   # standard EMBL column for key
            if key in ft_keep:
                current_ft = key
                keep_lines.append(l)
            if current_ft in ft_keep:
                subkey = l[21:].split('=')[0].strip('/')
                if subkey in ft_keep[current_ft]:
                    keep_lines.append(l)
        elif tag == 'SQ' and not sq_added:
            sq_added=True
            keep_lines.append(format_sq_block(str(consensus_seq)))
        else:   
            keep_lines.append(l)
    return "\n".join(keep_lines) + "\n"


def main():
    args = parse_args()

    # 1) Read all consensus sequences (segments) for this sample
    sys.stderr.write(f"[INFO] Reading consensus FASTA: {args.consensus}\n")
    with open(args.consensus) as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    if not records:
        sys.stderr.write("[ERROR] No sequences in consensus FASTA.\n")
        sys.exit(1)

    sys.stderr.write(f"[INFO] Found {len(records)} segment(s) in consensus.\n")

    # 2) Process each segment separately and collect patched EMBL blocks
    patched_records = []

    for idx, rec in enumerate(records, start=1):
        accession = rec.id.split()[0]
        sys.stderr.write(f"[INFO] Segment {idx}: accession={accession}\n")

        try:
            if accession.startswith("NC_"):
                embl_text = fetch_from_ncbi_genbank(accession)
            else:
                embl_text = fetch_embl_text(accession)
            patched = patch_single_embl(embl_text, args.sample_id, accession, rec.seq)
        except RuntimeError as e:
            sys.stderr.write(f"[ERROR] {e}\n")
            sys.exit(1)

        patched_records.append(patched.strip("\n"))

    # 3) Concatenate all EMBL records into one output file
    sys.stderr.write(f"[INFO] Writing {len(patched_records)} EMBL record(s) to {args.output}\n")
    with open(args.output, "w") as out_handle:
        out_handle.write("\n".join(patched_records) + "\n")

    sys.stderr.write("[INFO] Done.\n")


if __name__ == "__main__":
    main()

import os
import csv
import hashlib


def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, "").count(os.sep)
        indent = "│   " * (level - 1) if level > 0 else ""
        dir_basename = os.path.basename(root)
        if root != startpath:
            print("{}├── {}/".format(indent, dir_basename))
        subindent = "│   " * level
        file_indent = "{}├── ".format(subindent)
        for i, f in enumerate(files):
            if i == len(files) - 1 and not dirs:
                file_indent = "{}└── ".format(subindent)
            print("{}{}".format(file_indent, f))


def txt_to_list(input_file):
    with open(input_file, "r", encoding="utf-8-sig") as f:
        return f.readlines()


def read_csv_file(input_file):
    with open(input_file, "r", encoding="utf-8-sig") as f:
        dialect = csv.Sniffer().sniff(f.readline())
        f.seek(0)
        return list(csv.DictReader(f, dialect=dialect))


# https://www.quickprogrammingtips.com/python/how-to-calculate-sha256-hash-of-a-file-in-python.html
def compute_checksum(file_path):
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)

    return sha256_hash.hexdigest()

import os, sys
import shutil


def _module_reader(f):
    for line in f:
        line = line.strip()
        if line == ">>END_MODULE":
            return
        yield line


def parse_fastqc_data_modules(f):
    for line in f:
        line = line.strip()
        if line.startswith(">>") and line != ">>END_MODULE":
            yield line[2:], _module_reader(f)


def base_first_value(base):
    if "-" in base:
        return int(base.split("-")[0])
    else:
        return int(base)


def base_last_value(base):
    if "-" in base:
        return int(base.split("-")[1])
    else:
        return int(base)


def select_last_bad_base(module, threshold):
    bad_sequence_started = False
    bad_sequence_ended = False
    for line in module:
        if bad_sequence_ended:
            break
        if line.startswith("#"):
            continue

        line_parts = line.split("\t")
        base = line_parts[0]
        concentrations = tuple(float(x) for x in line_parts[1:])
        diff = max(concentrations) - min(concentrations)
        if diff > threshold:
            bad_sequence_started = True
        else:
            if bad_sequence_started:
                return base_first_value(base) - 1

    if not bad_sequence_started:
        return 0
    elif not bad_sequence_ended:
        return base_last_value(base)

def select_last_bad_base_for_txt(txt_file, threshold=6):
    with open(txt_file, 'r') as f:
        for module_name, module in parse_fastqc_data_modules(f):
            if module_name.startswith("Per base sequence content"):
                return select_last_bad_base(module, threshold)

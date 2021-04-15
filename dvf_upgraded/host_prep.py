from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import os
import math
from pathlib import Path
import fire


def host_prep(
        random_seed,
        fragment_length,
        saving_path,
        plant,
        plant_path,
        bact_path,
):
    random.seed(a=random_seed)
    plant_seq_dict = {file.strip(".fna"): [Path(plant_path, file), 0, 0] for file in os.listdir(plant_path) if
                      file.startswith("chr")}
    bact_seq_dict = {file.strip(".fna"): [Path(bact_path, file), 0, 0] for file in os.listdir(bact_path) if
                     file.endswith("fna")}

    # taking host val and train size from virus seq size.
    vir_train_size = Path(saving_path, f'train/vir_tr_{plant}_rs{random_seed}.fasta').stat().st_size
    vir_val_size = Path(saving_path, f'val/vir_val_{plant}_rs{random_seed}.fasta').stat().st_size
    # fragment_length = 1000
    host_train_count = int(vir_train_size / fragment_length)
    host_val_count = int(vir_val_size / fragment_length)
    # we select 70% of contigs from plants and 30% of contigs from bacteria
    plant_train_count = int(host_train_count * 0.7)
    plant_val_count = int(host_val_count * 0.7)
    bact_train_count = int(host_train_count * 0.3)
    bact_val_count = int(host_val_count * 0.3)

    for i in range(plant_train_count):
        seq = random.choice(list(plant_seq_dict.keys()))
        plant_seq_dict[seq][1] += 1
    for i in range(plant_val_count):
        seq = random.choice(list(plant_seq_dict.keys()))
        plant_seq_dict[seq][2] += 1
    for i in range(bact_train_count):
        seq = random.choice(list(bact_seq_dict.keys()))
        bact_seq_dict[seq][1] += 1
    for i in range(bact_val_count):
        seq = random.choice(list(bact_seq_dict.keys()))
        bact_seq_dict[seq][2] += 1

    # merging two host dictionaries
    seq_dict = {**plant_seq_dict, **bact_seq_dict}

    tr = []
    for seq_name in list(seq_dict.keys()):
        i = 0
        fasta = list(SeqIO.parse(seq_dict[seq_name][0], "fasta"))
        while i < seq_dict[seq_name][1]:
            # select chromosomes if there are any
            fasta_entry = random.choice(fasta)
            rand_end = len(fasta_entry.seq) - fragment_length
            start = random.randrange(rand_end)
            fragment = fasta_entry.seq[start:start + fragment_length]
            assert len(fragment) == fragment_length
            # fight gap bias
            if fragment.count('N') + fragment.count('-') < fragment_length * 0.05:
                record = SeqRecord(
                    seq=fragment,
                    id=f"{seq_name}_{fragment_length}_{start}",
                )
                tr.append(record)
                i += 1

    random.shuffle(tr)
    print(len(tr))
    SeqIO.write(tr, Path(saving_path, f"train/host_tr_{plant}_rs{random_seed}.fasta"), "fasta")

    val = []
    for seq_name in list(seq_dict.keys()):
        i = 0
        fasta = list(SeqIO.parse(seq_dict[seq_name][0], "fasta"))
        while i < seq_dict[seq_name][2]:
            # select chromosomes if there are any
            fasta_entry = random.choice(fasta)
            rand_end = len(fasta_entry.seq) - fragment_length
            start = random.randrange(rand_end)
            fragment = fasta_entry.seq[start:start + fragment_length]
            assert len(fragment) == fragment_length
            # fight gap bias
            if fragment.count('N') + fragment.count('-') < fragment_length * 0.05:
                record = SeqRecord(
                    seq=fasta_entry.seq[start:start + fragment_length],
                    id=f"{seq_name}_{fragment_length}_{start}",
                )
                val.append(record)
                i += 1

    random.shuffle(val)
    print(len(val))
    SeqIO.write(val, Path(saving_path, f"val/host_val_{plant}_rs{random_seed}.fasta"), "fasta")
    print("host preparation finished")

if __name__ == '__main__':
    fire.Fire(host_prep)


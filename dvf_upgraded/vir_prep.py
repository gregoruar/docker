from Bio import SeqIO
import random
import math
from pathlib import Path
import fire


def vir_prep(
        virus_path,
        random_seed,
        saving_path,
        plant,
):
    random.seed(a=random_seed)
    seqs_path = Path(virus_path)
    seqs = list(SeqIO.parse(seqs_path, "fasta"))
    random.shuffle(seqs)

    tr = []
    val = []
    threshold_80 = math.floor(len(seqs) * 0.8)
    for i in range(len(seqs)):
        if i < threshold_80:
            tr.append(seqs[i])
        else:
            val.append(seqs[i])
    print(len(tr))
    print(len(val))
    SeqIO.write(tr, Path(saving_path, f"train/vir_tr_{plant}_rs{random_seed}.fasta"), "fasta")
    SeqIO.write(val, Path(saving_path, f"val/vir_val_{plant}_rs{random_seed}.fasta"), "fasta")
    print("finished vir preparation")


if __name__ == '__main__':
    fire.Fire(vir_prep)

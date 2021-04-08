from Bio import SeqIO
import random
import math
from pathlib import Path
import fire


def vir_prep(
        random_seed,
        saving_path,
):
    random.seed(a=random_seed)
    family = Path(saving_path).name
    seqs_path = Path(
        f"/mnt/cbib/inextvir/workspace/gregoruar/deepvirfinder/dvf_train/vir_db/families/plant-vir_NO-{family}_2020-01-04.fasta")
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
    SeqIO.write(tr, Path(saving_path, f"train/vir_tr_rs{random_seed}.fasta"), "fasta")
    SeqIO.write(val, Path(saving_path, f"val/vir_val_rs{random_seed}.fasta"), "fasta")


if __name__ == '__main__':
    fire.Fire(vir_prep)


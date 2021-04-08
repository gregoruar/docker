import numpy as np
from pathlib import Path
from sklearn.utils import shuffle
import re
import fire


def undersampling(
        random_seed,
        saving_path,
):
    for ds in ["train", "val"]:
        path = Path(saving_path, f"{ds}/encode")
        path.rename(path.parent / f"encode_rs{random_seed}")
        if ds == "train":
            ds = "tr"
        file_list = path.glob(f"host#host_{ds}#1k_num1_seq*")
        for p_ in file_list:
            print(Path(p_.parent, f"{p_.stem}_old{p_.suffix}"))
            p_.rename(Path(p_.parent, f"{p_.stem}_old{p_.suffix}"))
        vir_file = list(path.glob(f"virus#vir_{ds}#1k_num1_seq*_codebw.npy"))[0]
        n_seq = int(re.search(rf"(?<=seq)[0-9]+", str(vir_file)).group(0))
        host = list(path.glob(f"host#host_{ds}#1k_num1_seq*_old.npy"))
        host_bw = np.load(host[0])
        host_fw = np.load(host[1])
        host_bw_n, host_fw_n = shuffle(host_bw, host_fw, random_state=random_seed, n_samples=n_seq)
        print(f"host#host_{ds}#1k_num1_seq{n_seq}_codebw.npy")
        np.save(Path(path, f"host#host_{ds}#1k_num1_seq{n_seq}_codebw.npy"), host_bw_n)
        np.save(Path(path, f"host#host_{ds}#1k_num1_seq{n_seq}_codefw.npy"), host_fw_n)
    # families = [
    #     "Caulimoviridae",
    #     "Closteroviridae",
    #     "Phenuiviridae",
    #     "Unclassified",
    # ]
    # path_encode = Path(saving_path, "encode")
    # path_encode.rename(path_encode.parent / f"encode_rs{random_seed}")
    # renaming files
    # for f in families:
    # # for fam in families:
    # for ds in ["train", "val"]:
    #     path = Path(saving_path, f"{ds}/encode_rs{random_seed}")


if __name__ == '__main__':
    fire.Fire(undersampling)


import numpy as np
from pathlib import Path
from sklearn.utils import shuffle
import re
import fire


def undersampling(random_seed, saving_path, plant,):
    for ds in ["train", "val"]:
        path_old = Path(saving_path, f"{ds}/encode")
        path = Path(saving_path, f"{ds}/encode_{plant}_rs{random_seed}")
        path.mkdir(exist_ok=True)
        file_list = path_old.glob(f"*{plant}_rs{random_seed}*")
        for p_ in file_list:
            p_.replace(path / p_.name)
        if ds == "train":
            ds = "tr"
        if not list(path.glob(f"host#host_{ds}_{plant}_rs{random_seed}#1k_num1_seq*_old.npy")):
            file_list = path.glob(f"host#host_{ds}_{plant}_rs{random_seed}#1k_num1_seq*")
            for p_ in file_list:
                p_.rename(Path(p_.parent, f"{p_.stem}_old{p_.suffix}"))
        vir_file = list(path.glob(f"virus#vir_{ds}_{plant}_rs{random_seed}#1k_num1_seq*_codebw.npy"))[0]
        n_seq = int(re.search(rf"(?<=seq)[0-9]+", str(vir_file)).group(0))
        host = list(path.glob(f"host#host_{ds}_{plant}_rs{random_seed}#1k_num1_seq*_old.npy"))
        host_bw = np.load(host[0])
        host_fw = np.load(host[1])
        host_bw_n, host_fw_n = shuffle(host_bw, host_fw, random_state=random_seed, n_samples=n_seq)
        np.save(Path(path, f"host#host_{ds}_{plant}_rs{random_seed}#1k_num1_seq{n_seq}_codebw.npy"), host_bw_n)
        np.save(Path(path, f"host#host_{ds}_{plant}_rs{random_seed}#1k_num1_seq{n_seq}_codefw.npy"), host_fw_n)
    print("finished undersampling")


if __name__ == '__main__':
    fire.Fire(undersampling)

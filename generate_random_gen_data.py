#!/usr/bin/env python3

import sys
import random

def write_random_gen_data_file(num_snps, num_samples, fn_out):
    snp_ids = random.sample(range(0, num_snps * 10), num_samples)
    ALLELES = ["A", "T", "C", "G"]
    with open(fn_out, "w") as f_out:
        for snp_id in snp_ids:
            a0 = ALLELES[random.randint(0, 3)]
            a1 = list(set(ALLELES).difference([a0]))[random.randint(0, 2)]
            row = [snp_id, f"rs{snp_id}", snp_id * 100, a0, a1]
            for j in range(0, num_samples):
                p0 = random.uniform(0.85, 1)
                p1 = random.uniform(0, 1 - p0)
                p2 = 1 - p0 - p1
                row.extend(map(lambda _: round(_, 2), random.sample([p0, p1, p2], 3)))

            print(" ".join(map(str, row)), file = f_out)


def write_random_sample_file(num_samples, fn_out):
    with open(fn_out, "w") as f_out:
        samples = random.sample(range(0, num_samples * 10), num_samples)
        max_sample_id_len = len(str(max(samples)))
        samples = ["sample_{}".format(str(_).zfill(max_sample_id_len)) for _ in samples]
        print("ID_1 ID_2 missing father mother sex plink_pheno", file = f_out)
        print("0 0 0 D D D B", file = f_out)
        print("\n".join([f"0 {_} 0 0 0 0 -9" for _ in samples]), file = f_out)


def main():
    if len(sys.argv) != 4:
        print(f"USAGE: {sys.argv[0]} num_snps, num_samples, prefix", file = sys.stderr)
        sys.exit(1)

    num_snps = int(sys.argv[1])
    num_samples = int(sys.argv[2])
    prefix = sys.argv[3]

    write_random_gen_data_file(num_snps, num_samples, f"{prefix}.gen")
    write_random_sample_file(num_samples, f"{prefix}.samples")


if __name__ == "__main__":
    main()

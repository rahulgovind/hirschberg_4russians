import sys
import argparse
import random
from datetime import datetime

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate random sequences')
    parser.add_argument("-n", dest='N', default=100, type=int,
                        action="store", help="Size of sequence")
    parser.add_argument("-c", dest="C", default=1, type=int,
                        action="store", help="Number of sequences to generate")
    parser.add_argument("-q", dest='Q', default="ATCG",
                        action="store", help="Character set")
    parser.add_argument("--seed", action="store", type=int, dest="seed",
                        help="Seed for randomization")
    parser.add_argument("-o", action="store", dest="outfile",
                        help="Output file to store")
    args = parser.parse_args(sys.argv[1:])

    seed = random.seed(datetime.now()) if args.seed is None else args.seed
    random.seed(seed)
    s = []
    for i in range(args.C):
        seq = [random.choice(args.Q) for _ in range(args.N)]
        s.append(''.join(seq))

    if args.outfile is not None:
        with open(args.outfile, 'w') as f:
            f.write("\n".join(s) + "\n")
    else:
        for si in s:
            print(si)

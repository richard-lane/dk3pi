import numpy.random as rnd
import sys
import os


def main(arg_list):
    with open("event_nums.txt", "w") as f:
        f.write(str(rnd.poisson(int(os.environ["MEANNUMDCS"]))) + "\n")
        f.write(str(rnd.poisson(int(os.environ["MEANNUMCF"]))))


if __name__ == "__main__":
    main(sys.argv)

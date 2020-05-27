import numpy as np

def main():
    r_pulls = []
    with open("pull.txt", "r") as f:
        fit_results = [line.rstrip() for line in f]

        for result in fit_results:
            result = result.split(',')
            r = result[0]
            dr = result[1]
            r_pulls.append((float(r) - 0.0549)/float(dr))
    print(np.mean(r_pulls), "+-", np.std(r_pulls))

if __name__ == '__main__':
    main()

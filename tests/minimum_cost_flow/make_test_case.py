from pathlib import Path


def aoj_grl_6_b(input_directory_path: Path):
    for input_file_path in input_directory_path.iterdir():
        if "in.in" not in str(input_file_path):
            continue

        expected_file_path = input_file_path.with_suffix(".out")
        with expected_file_path.open("r") as f:
            expected = int(f.readline().strip())
            if expected == -1:
                expected = "infeasible"

        lines = []
        with input_file_path.open("r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    n, m, f = map(int, line.strip().split())
                    lines.append(f"{n} {m} {expected}")
                    b = [0] * n
                    (b[0], b[-1]) = (f, -f)
                    lines.append("\n".join([str(x) for x in b]))
                else:
                    f, t, u, c = line.strip().split()
                    lines.append(f"{f} {t} {0} {u} {c}")

        output_file_path = Path(str(input_file_path).replace(".in.in", ".txt"))
        with output_file_path.open("w") as f:
            f.write("\n".join(lines))

        print(input_file_path, expected_file_path)


def library_checker_min_cost_b_flow(input_directory_path: Path):
    for input_file_path in input_directory_path.iterdir():
        if ".in" not in str(input_file_path):
            continue

        expected_file_path = input_file_path.with_suffix(".out")
        with expected_file_path.open("r") as f:
            expected = f.readline().strip()
            if expected == -1:
                expected = "infeasible"

        lines = []
        with input_file_path.open("r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    n, m = map(int, line.strip().split())
                    lines.append(f"{n} {m} {expected}")
                else:
                    lines.append(line.strip())

        output_file_path = Path(str(input_file_path).replace(".in", ".txt"))
        with output_file_path.open("w") as f:
            f.write("\n".join(lines))


def main():
    aoj_grl_6_b(Path("AOJ_GRL_6_B"))
    library_checker_min_cost_b_flow(Path("LibraryChecker_min_cost_b_flow"))


if __name__ == "__main__":
    main()

from pathlib import Path


def aoj_grl_6_b(input_directory_path: Path):
    for input_file_path in input_directory_path.iterdir():
        if "in.in" not in str(input_file_path):
            continue

        expected_file_path = input_file_path.with_suffix(".out")
        with expected_file_path.open("r") as f:
            expected = int(f.readline().strip())

        lines = []
        with input_file_path.open("r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    n, m = map(int, line.strip().split())
                    source, sink = 0, n - 1
                    lines.append(f"{n} {m} {source} {sink} {expected}")
                else:
                    f, t, u = line.strip().split()
                    lines.append(f"{f} {t} {u}")

        output_file_path = Path(str(input_file_path).replace(".in.in", ".txt"))
        with output_file_path.open("w") as f:
            f.write("\n".join(lines))


def libre_oj_101(input_directory_path: Path):
    for input_file_path in input_directory_path.iterdir():
        if ".in" not in str(input_file_path):
            continue

        expected_file_path = input_file_path.with_suffix(".out")
        with expected_file_path.open("r") as f:
            expected = int(f.readline().strip())

        lines = []
        with input_file_path.open("r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    n, m, source, sink = map(int, line.strip().split())
                    source -= 1
                    sink -= 1
                    lines.append(f"{n} {m} {source} {sink} {expected}")
                else:
                    f, t, u = map(int, line.strip().split())
                    lines.append(f"{f - 1} {t - 1} {u}")

        output_file_path = Path(str(input_file_path).replace(".in", ".txt"))
        with output_file_path.open("w") as f:
            f.write("\n".join(lines))


def main():
    aoj_grl_6_b(Path("AOJ_GRL_6_A"))
    libre_oj_101(Path("LibreOJ_101"))


if __name__ == "__main__":
    main()

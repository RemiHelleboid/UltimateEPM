from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
import os
import vtk


def convert_vtp(input_path: Path, output_path: Path, appended: bool = True) -> None:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(input_path))
    reader.Update()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(output_path))
    writer.SetInputData(reader.GetOutput())

    if appended:
        writer.SetDataModeToAppended()
        writer.EncodeAppendedDataOff()
    else:
        writer.SetDataModeToBinary()

    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write: {output_path}")


def convert_one(task: tuple[str, str, bool]) -> tuple[str, str]:
    input_path = Path(task[0])
    output_path = Path(task[1])
    appended = task[2]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    convert_vtp(input_path=input_path, output_path=output_path, appended=appended)
    return str(input_path), str(output_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_folder")
    parser.add_argument("output_folder")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--inline-binary", action="store_true")
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel worker processes. Use 0 to use all available CPUs.",
    )
    args = parser.parse_args()

    input_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)
    jobs = (os.cpu_count() or 1) if args.jobs == 0 else args.jobs
    if jobs < 1:
        raise ValueError("--jobs must be positive, or 0 to use all available CPUs.")

    pattern = "**/*.vtp" if args.recursive else "*.vtp"
    files = sorted(input_folder.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No .vtp files found in {input_folder}")

    tasks: list[tuple[str, str, bool]] = []
    for input_path in files:
        relative_path = input_path.relative_to(input_folder)
        output_path = output_folder / relative_path
        tasks.append((str(input_path), str(output_path), not args.inline_binary))

    if jobs == 1:
        results = [convert_one(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            results = list(executor.map(convert_one, tasks))

    for input_path, output_path in results:
        print(f"Converted: {input_path} -> {output_path}")

    print(f"Done. Converted {len(files)} file(s).")


if __name__ == "__main__":
    main()

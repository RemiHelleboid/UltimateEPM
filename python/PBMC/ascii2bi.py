from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
import os

import vtk


SUPPORTED_EXTENSIONS = {".vtp", ".vtu"}


def create_reader(input_path: Path):
    extension = input_path.suffix.lower()

    if extension == ".vtp":
        return vtk.vtkXMLPolyDataReader()

    if extension == ".vtu":
        return vtk.vtkXMLUnstructuredGridReader()

    raise ValueError(f"Unsupported VTK file extension: {input_path.suffix}")


def create_writer(output_path: Path):
    extension = output_path.suffix.lower()

    if extension == ".vtp":
        return vtk.vtkXMLPolyDataWriter()

    if extension == ".vtu":
        return vtk.vtkXMLUnstructuredGridWriter()

    raise ValueError(f"Unsupported VTK file extension: {output_path.suffix}")


def convert_vtk_file(
    input_path: Path,
    output_path: Path,
    appended: bool = True,
) -> None:
    reader = create_reader(input_path)
    reader.SetFileName(str(input_path))
    reader.Update()

    output = reader.GetOutput()

    if output is None:
        raise RuntimeError(f"Failed to read: {input_path}")

    writer = create_writer(output_path)
    writer.SetFileName(str(output_path))
    writer.SetInputData(output)

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

    convert_vtk_file(
        input_path=input_path,
        output_path=output_path,
        appended=appended,
    )

    return str(input_path), str(output_path)


def discover_files(input_folder: Path, recursive: bool) -> list[Path]:
    iterator = input_folder.rglob("*") if recursive else input_folder.glob("*")

    return sorted(
        path
        for path in iterator
        if path.is_file() and path.suffix.lower() in SUPPORTED_EXTENSIONS
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert VTP and VTU files to appended or inline binary XML."
    )
    parser.add_argument("input_folder", type=Path)
    parser.add_argument("output_folder", type=Path)
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--inline-binary", action="store_true")
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Worker processes. Use 0 for all available CPUs.",
    )
    args = parser.parse_args()

    input_folder = args.input_folder.resolve()
    output_folder = args.output_folder.resolve()

    if not input_folder.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_folder}")

    if input_folder == output_folder:
        raise ValueError("Input and output directories must be different.")

    jobs = (os.cpu_count() or 1) if args.jobs == 0 else args.jobs

    if jobs < 1:
        raise ValueError("--jobs must be positive, or 0 to use all CPUs.")

    files = discover_files(input_folder, args.recursive)

    if not files:
        supported = ", ".join(sorted(SUPPORTED_EXTENSIONS))
        raise FileNotFoundError(
            f"No supported VTK files ({supported}) found in {input_folder}"
        )

    tasks = [
        (
            str(input_path),
            str(output_folder / input_path.relative_to(input_folder)),
            not args.inline_binary,
        )
        for input_path in files
    ]

    if jobs == 1:
        results = [convert_one(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            results = list(executor.map(convert_one, tasks))

    for input_path, output_path in results:
        print(f"Converted: {input_path} -> {output_path}")

    print(f"Done. Converted {len(results)} file(s).")


if __name__ == "__main__":
    main()
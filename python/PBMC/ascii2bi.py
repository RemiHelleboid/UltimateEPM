from pathlib import Path
import argparse
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_folder")
    parser.add_argument("output_folder")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--inline-binary", action="store_true")
    args = parser.parse_args()

    input_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)

    pattern = "**/*.vtp" if args.recursive else "*.vtp"
    files = sorted(input_folder.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No .vtp files found in {input_folder}")

    for input_path in files:
        relative_path = input_path.relative_to(input_folder)
        output_path = output_folder / relative_path
        output_path.parent.mkdir(parents=True, exist_ok=True)

        convert_vtp(
            input_path=input_path,
            output_path=output_path,
            appended=not args.inline_binary,
        )

        print(f"Converted: {input_path} -> {output_path}")

    print(f"Done. Converted {len(files)} file(s).")


if __name__ == "__main__":
    main()
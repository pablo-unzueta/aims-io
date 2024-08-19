import pandas as pd
import csv


def create_csv_from_aims_traj_dump(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)

        # Read the header line, remove the "#" character, and split it into columns
        header = infile.readline().strip().lstrip("#").split()
        writer.writerow(header)

        # Read the data lines and split each into columns
        for line in infile:
            # Split the line into individual components
            data = line.strip().split()
            # Convert all string data to float or int as needed
            data = [float(item) if "." in item else int(item) for item in data]
            writer.writerow(data)


def test_create_csv_from_aims_traj_dump():
    create_csv_from_aims_traj_dump(
        "/Users/pablo/test-aims/0000/TrajDump.1", "/Users/pablo/test-aims/0000/TrajDump.1.csv"
    )


if __name__ == "__main__":
    test_create_csv_from_aims_traj_dump()

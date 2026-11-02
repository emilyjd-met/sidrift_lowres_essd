import os
import argparse
from argparse import RawDescriptionHelpFormatter
from datetime import datetime, timedelta
import netCDF4 as nc
import shutil  # For copying files


def parse_args():
    p = argparse.ArgumentParser("Divide matchup files into two sets based on time",
                                 formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--indir', required=True,
                   help="Input directory containing the matchup files")
    p.add_argument('-a', '--outdir_a', required=True,
                   help="Output directory A (set A)")
    p.add_argument('-b', '--outdir_b', required=True,
                   help="Output directory B (set B)")
    p.add_argument('-s', '--start_date', required=True,
                   help="Start date in YYYYMMDD format")
    p.add_argument('-e', '--end_date', required=True,
                   help="End date in YYYYMMDD format")
    p.add_argument('-z', '--hemi', required=True,
                   help="Hemisphere")
    p.add_argument('-t', '--inst', required=True,
                   help="Instrument name")
    p.add_argument('-c', '--channel', required=True,
                   help="Channel name")
    p.add_argument('-d', '--duration', type=int, required=True,
                   help="Duration in days before the processing date")
    p.add_argument('-q', '--split_duration', type=int, default=7,
                   help="Number of days for each setA and setB cycle (default: 7 days)")

    args = p.parse_args()
    return args


def get_matchup_size(filepath):
    """
    Reads the size of the 'matchup' dimension from the NetCDF file.
    """
    try:
        with nc.Dataset(filepath, 'r') as ds:
            if 'matchup' in ds.dimensions:
                return len(ds.dimensions['matchup'])
            else:
                print(f"Warning: 'matchup' dimension not found in {filepath}")
                return 0
    except Exception as e:
        print(f"Failed to read NetCDF file {filepath}: {e}")
        return 0


def copy_file(filepath, output_dir, date):
    """
    Copies the file to the specified output directory, maintaining the
    YYYY/MM subdirectory structure.
    """
    # Create the YYYY/MM subdirectory based on the provided date
    yyyy_mm_subdir = os.path.join(output_dir, date.strftime("%Y/%m"))
    os.makedirs(yyyy_mm_subdir, exist_ok=True)

    # Copy the file to the appropriate subdirectory
    destination = os.path.join(yyyy_mm_subdir, os.path.basename(filepath))
    shutil.copy(filepath, destination)
    print(f"Copied {filepath} to {destination}")


def split_matchups(indir, outdir_a, outdir_b, start_date, end_date, hemi,
                   inst, channel, duration, split_duration):
    """
    Cycles over each date between start_date and end_date, divides files
    into setA or setB based on the YYYYMMDD in their file pattern, and
    processes the files.
    """
    # Convert start_date and end_date to datetime objects
    start_date = datetime.strptime(start_date, "%Y%m%d")
    end_date = datetime.strptime(end_date, "%Y%m%d")

    # Initialize counters
    setA_count = 0
    setB_count = 0
    setA_matchup_total = 0
    setB_matchup_total = 0

    # Cycle through each date in the range
    current_date = start_date
    while current_date <= end_date:
        # Calculate the YYYYMMDD for duration days before the current
        # processing date
        x_date = current_date - timedelta(days=duration)
        current_yyyymmdd = current_date.strftime("%Y%m%d")
        x_yyyymmdd = x_date.strftime("%Y%m%d")

        # Construct the file pattern
        file_pattern = f"{indir}/{x_date.year:04d}/{x_date.month:02d}/matchup_{inst}_{channel}_simplex_lev2_{hemi}-ease2-750_{x_yyyymmdd}12w_{current_yyyymmdd}12w_NN3Dcol.nc"

        # Check if the file exists
        if os.path.isfile(file_pattern):
            # Determine whether the file belongs to set A or set B based on
            # the second date (YYYYMMDD)
            file_date = current_yyyymmdd
            file_datetime = datetime.strptime(file_date, "%Y%m%d")
            days_since_start = (file_datetime - start_date).days
            # Alternate every split_duration days between set A and set B
            if (days_since_start // split_duration) % 2 == 0:
                output_dir = outdir_a
                setA_count += 1
            else:
                output_dir = outdir_b
                setB_count += 1

            # Copy the file to the appropriate directory
            copy_file(file_pattern, output_dir, current_date)

            # Get the size of the matchup dimension and add it to the
            # appropriate total
            matchup_size = get_matchup_size(file_pattern)
            if output_dir == outdir_a:
                setA_matchup_total += matchup_size
            else:
                setB_matchup_total += matchup_size
        else:
            print(f"File not found: {file_pattern}")
        
        # Move to the next date
        current_date += timedelta(days=1)

    # Print the results
    print("\n=== Processing Results ===")
    print(f"Set A: {setA_count} files, {setA_matchup_total} matchups")
    print(f"Set B: {setB_count} files, {setB_matchup_total} matchups")


def main():
    args = parse_args()

    # Call the main processing function
    split_matchups(
        indir=args.indir,
        outdir_a=args.outdir_a,
        outdir_b=args.outdir_b,
        start_date=args.start_date,
        end_date=args.end_date,
        hemi=args.hemi,
        inst=args.inst,
        channel=args.channel,
        duration=args.duration,
        split_duration=args.split_duration
    )


if __name__ == '__main__':
    main()

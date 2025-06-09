import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import os
import sys


input_bdg = sys.argv[1]

# Read data
bedgraph = pd.read_csv(input_bdg, sep="\t", header=None, names=["chrom", "start", "end", "value"], skiprows=1)
#bedgraph = bedgraph[bedgraph['chrom'] == 'chr19']
finer_resolution_step = 10
# Output file
output_file = sys.argv[2]

# Write the track line to the output file
with open(output_file, "w") as f:
    f.write("track type=bedGraph\n")

    # Only process chromosome 19
    if not bedgraph.empty:
        # Initialize the start index
        i = 0
        
        while i < len(bedgraph):
            row = bedgraph.iloc[i]
            
            if row['value'] == 0.0:
                # If the value is zero, write the row as-is
                f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['value']}\n")
                i += 1
            else:
                # Find the start of the non-zero segment
                start = i
                # Move to the end of the non-zero segment
                while i < len(bedgraph) and bedgraph.iloc[i]['value'] != 0.0:
                    i += 1
                end = i

                # Get the segment data
                non_zero_segment = bedgraph.iloc[start:end]

                # Use only start positions for x, and corresponding y values
                x = non_zero_segment['start'].values
                y = non_zero_segment['value'].values

                # Check if x and y have at least two points
                if len(x) > 1:
                    # Interpolate at a finer resolution for the non-zero segment
                    interpolated_x = np.arange(x.min(), x.max(), finer_resolution_step)
                    f_interp = interp1d(x, y, kind='linear', fill_value="extrapolate")
                    interpolated_y = f_interp(interpolated_x)

                    # Combine to a DataFrame for this non-zero segment
                    interpolated_frame = pd.DataFrame({
                        "chrom": non_zero_segment.iloc[0]['chrom'],
                        "start": interpolated_x[:-1],
                        "end": interpolated_x[1:],  # Ensure no overlapping
                        "value": interpolated_y[:-1]
                    })

                    # Append the interpolated data to the file
                    interpolated_frame.to_csv(f, sep="\t", header=False, index=False, mode='a')
                else:
                    # If not enough data points, write the row directly
                    for j in range(len(non_zero_segment)):
                        f.write(f"{non_zero_segment.iloc[j]['chrom']}\t{non_zero_segment.iloc[j]['start']}\t{non_zero_segment.iloc[j]['end']}\t{non_zero_segment.iloc[j]['value']}\n")
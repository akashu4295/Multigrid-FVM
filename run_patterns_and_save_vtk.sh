#!/bin/bash

echo "Compiling the code on GPU..."
if nvc -acc -Minfo=accel Version-6-26-04.c; then
    echo "GPU compilation successful"
else
    echo "GPU compilation failed"
    echo "Compiling for CPU instead..."
    nvc Version-6-26-04.c 
fi

reynolds_numbers=("20")
patterns_dir="patterns"

# Grid size configuration for 16x16 grid with dynamic column calculations
base_n=16
base_col1=$(( 5 * base_n * (2 ** (6 - 1)) ))
base_col2=$(( base_n * (2 ** (6 - 1)) ))
base_col3=6
base_col4=6
base_col5="$base_n"

for reynolds in "${reynolds_numbers[@]}"; do
    inv_re=$(echo "scale=6; 1 / $reynolds" | bc)
    echo "Setting up InputData for Reynolds Number: Re = $reynolds (1/Re = $inv_re)"
    
    # Update inverse Reynolds number in Row 2, Column 3
    awk -v new_val="$inv_re" 'NR==2 {$3=new_val} {print}' InputData > InputData.tmp && mv InputData.tmp InputData

    # Update Line 1 of InputData with calculated dimensions
    awk -v c1="$base_col1" -v c2="$base_col2" -v c3="$base_col3" -v c4="$base_col4" -v c5="$base_col5" \
        'NR==1 {$1=c1; $2=c2; $3=c3; $4=c4; $5=c5} {print}' InputData > InputData.tmp && mv InputData.tmp InputData
    
    echo "--> Configured InputData for N = $base_n -> Row 1: [$base_col1 $base_col2 $base_col3 $base_col4 $base_col5]"

    # Build destination directories under solution/Re_<reynolds>/
    sol_vtk_dir="solution/Re_${reynolds}/vtk_files"
    sol_conv_dir="solution/Re_${reynolds}/convergence_files"
    mkdir -p "$sol_vtk_dir"
    mkdir -p "$sol_conv_dir"

    # Process all pattern_*.txt files directly inside the patterns folder
    for filepath in "${patterns_dir}"/pattern_*.txt; do
        
        [ -e "$filepath" ] || continue

        filename=$(basename "$filepath")
        number="${filename#pattern_}"
        number="${number%.txt}"
        
        echo "  -> Running $filename for Re = $reynolds..."

        # Copy active pattern to execution directory as pattern.txt
        cp "$filepath" ./pattern.txt

        # Run simulation
        ./a.out

        # Save convergence output
        if [ -f "convergence.csv" ]; then
            mv convergence.csv "${sol_conv_dir}/convergence_${number}.csv"
        else
            echo "  [Error] ./a.out did not produce convergence.csv for $filename"
        fi

        # Save VTK output
        if [ -f "vtk_acc.vtk" ]; then
            cp vtk_acc.vtk "${sol_vtk_dir}/vtk_${number}.vtk"
            rm vtk_acc.vtk
        else
            echo "  [Error] ./a.out did not produce vtk_acc.vtk for $filename"
        fi

        # Clean up temporary pattern file
        rm -f ./pattern.txt

    done
done

echo "--------------------------------------------"
echo "All loops completed. Running post-processing..."
echo "Removing files that haven't converged...."

python remove_files.py

echo "All processing completed successfully."
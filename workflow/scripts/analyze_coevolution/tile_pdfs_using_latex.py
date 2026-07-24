#!/usr/bin/env python3
"""
Generate LaTeX files using TikZ to tile coevolution analysis PDFs.
Creates a single page with 8 rows x 3 columns layout.
"""

import os
import subprocess

# Define the MAFL combinations in order
MAFL_COMBINATIONS = [
    "activation",
    "maf_tf",
    "ehr_basic_bzip",
    "basic_bzip",
    "bzip",
    "ehr",
    "basic",
    "ehr_basic",
]

# Add labels for rows (optional)
row_labels = {
    "activation": "Activation",
    "maf_tf": "MAF TF",
    "ehr_basic_bzip": "EHR Basic bZIP",
    "basic_bzip": "Basic bZIP",
    "bzip": "bZIP",
    "ehr": "EHR",
    "basic": "Basic",
    "ehr_basic": "EHR Basic",
}

# Define the species variants
SPECIES = ["all", "mammal", "nonmammal"]

# Base directory
BASE_DIR = snakemake.params["BASE_DIR"] #"scratch/coevolution_analysis"
print(BASE_DIR)

otx_domain = snakemake.wildcards["otxdom"]

def get_pdf_path(otx_type, mafl_combo, species):
    """Generate the PDF file path."""
    return f"coevol_summary__OTX__{otx_type}__MAFL__{mafl_combo}__{species}.pdf"

def create_latex_document(otx_type, output_tex):
    """
    Create a LaTeX document using TikZ to tile PDFs.

    Args:
        otx_type: Either 'homeobox' or 'otx1_tf'
        output_tex: Path for the output .tex file
    """

    latex_content = r"""\documentclass[border=0pt]{standalone}
\usepackage{graphicx}
\usepackage{tikz}

\begin{document}
\begin{tikzpicture}

"""

    # Define spacing and scale
    # Adjust these values based on your PDF sizes
    col_width = 6.0  # cm between columns
    row_height = 4.0  # cm between rows
    scale = 0.38  # scale factor for PDFs

    print(f"Generating LaTeX for {otx_type}...")

    for row_idx, mafl_combo in enumerate(MAFL_COMBINATIONS):
        for col_idx, species in enumerate(SPECIES):
            pdf_file = get_pdf_path(otx_type, mafl_combo, species)
            pdf_path = f"../{pdf_file}"

            # Calculate position (origin at top-left)
            x_pos = col_idx * col_width
            y_pos = -row_idx * row_height

            # Check if file exists
            full_path = os.path.join(BASE_DIR, pdf_file)
            if not os.path.exists(full_path):
                print(f"  WARNING: Missing {pdf_file}")
                continue

            # Add node with PDF
            latex_content += f"  \\node[anchor=north west] at ({x_pos},{y_pos}) {{\n"
            latex_content += f"    \\includegraphics[scale={scale},page=1]{{{pdf_path}}}\n"
            latex_content += f"  }};\n"
            y_pos = -row_idx * row_height - row_height/2
            label = row_labels[mafl_combo]
            latex_content += f"  \\node[anchor=east] at (-0.3,{y_pos}) {{\\small {label}}};\n"

    # Add labels for columns (optional)
    latex_content += r"""
  % Column labels
  \node[anchor=south] at (0,0.3) {\small All};
  \node[anchor=south] at (6.5,0.3) {\small Mammal};
  \node[anchor=south] at (13,0.3) {\small Non-mammal};

"""

    latex_content += r"""
\end{tikzpicture}
\end{document}
"""

    # Write the LaTeX file
    with open(output_tex, 'w') as f:
        f.write(latex_content)

    print(f"Created: {output_tex}")
    return output_tex

def compile_latex(tex_file):
    """Compile LaTeX file to PDF using pdflatex."""
    print(f"\nCompiling {tex_file}...")

    # Get directory and filename
    tex_dir = os.path.dirname(tex_file)

    try:
        # Run pdflatex twice for proper positioning
        for i in range(2):
            result = subprocess.run(
                ['pdflatex', '-interaction=nonstopmode', os.path.basename(tex_file)],
                cwd=tex_dir,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                print(f"Error compiling LaTeX (run {i+1}):")
                print(result.stdout)
                return False

        pdf_file = tex_file.replace('.tex', '.pdf')
        print(f"Successfully created: {pdf_file}")

        # Clean up auxiliary files
        for ext in ['.aux', '.log']:
            aux_file = tex_file.replace('.tex', ext)
            if os.path.exists(aux_file):
                os.remove(aux_file)

        return True

    except FileNotFoundError:
        print("ERROR: pdflatex not found. Please install a LaTeX distribution (e.g., TeX Live, MiKTeX)")
        return False

def main():
    """Main function to create both tiled PDFs."""
    # Create output directory if it doesn't exist
    output_dir = f"{BASE_DIR}/tiled"
    os.makedirs(output_dir, exist_ok=True)

    # Create LaTeX and compile for homeobox
    print("=" * 60)
    print("Creating homeobox tiled PDF...")
    print("=" * 60)
    tex_file = create_latex_document(otx_domain, f"{output_dir}/OTX_{otx_domain}_MAFL_tiled.tex")
    compile_latex(tex_file)

    print("\n" + "=" * 60)
    print("Done! Tiled PDFs created in:", output_dir)
    print("=" * 60)
    print("\nYou can adjust the 'scale', 'col_width', and 'row_height'")
    print("variables in the script to fine-tune the layout.")

if __name__ == "__main__":
    main()

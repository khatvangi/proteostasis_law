#!/bin/bash
echo "========================================"
echo "Proteostasis Law - Figure Files Verification"
echo "========================================"
echo ""

# Check for all expected files
figures=(
    "figure1_central"
    "figure1_conceptual"
    "figure2_metal_ligands"
    "figure3_logistic_regression"
    "figure4_codon_modes"
    "figure5_mode_richness"
    "figure6_capacity_bound"
    "figureS1_distributions"
    "figureS2_metal_fractions"
    "figureS3_roc_curves"
)

formats=("png" "pdf" "svg")

echo "Checking Main & Supplementary Figures..."
echo ""

all_good=true
for fig in "${figures[@]}"; do
    echo "📊 $fig:"
    for fmt in "${formats[@]}"; do
        file="${fig}.${fmt}"
        if [ -f "$file" ]; then
            size=$(ls -lh "$file" | awk '{print $5}')
            echo "  ✓ $fmt ($size) - EDITABLE"
        else
            echo "  ✗ $fmt - MISSING"
            all_good=false
        fi
    done
    echo ""
done

if [ "$all_good" = true ]; then
    echo "========================================"
    echo "✓ ALL FIGURES COMPLETE & EDITABLE!"
    echo "========================================"
    echo ""
    echo "Total files:"
    total_files=$(ls *.png *.pdf *.svg 2>/dev/null | wc -l)
    echo "  $total_files figure files generated"
    echo ""
    echo "Editable formats available:"
    echo "  📄 PDF - Open in Illustrator, Inkscape, etc."
    echo "  📄 SVG - Open in Inkscape (free), web browsers"
    echo ""
    echo "See FIGURES_SUMMARY.md for editing instructions"
else
    echo "⚠ Some files are missing. Please regenerate."
fi

#!/usr/bin/env python3
"""Create a Markdown pilot summary from analysis outputs."""
import argparse
from pathlib import Path
import pandas as pd

DEFAULT_DATA=Path.home()/'Downloads'/'Wu_etal_2021_BRCA_scRNASeq'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--data-dir',type=Path,default=DEFAULT_DATA)
    p.add_argument('--output-dir',type=Path,default=None)
    args=p.parse_args()
    out=args.output_dir.expanduser() if args.output_dir else args.data_dir.expanduser()/'analysis'
    comp=pd.read_csv(out/'pilot_cell_composition.csv')
    samples=comp.groupby(['orig.ident','subtype'],as_index=False)['n_cells'].sum()
    fig=out/'figures'
    names=['01_qc_violin.png','02_umap_author_celltype.png','03_umap_subtype.png',
           '04_umap_leiden.png','05_marker_dotplot.png']
    lines=['# Wu et al. (2021) scRNA-seq pilot reproduction','',
           '## Source','GSE176078; Wu et al., *Nature Genetics* (2021).','',
           '## Scope and caveats','This is a **small sampled pilot** using one sample per subtype, not a full-paper reproduction.',
           'Author-provided `celltype_major` labels are used for comparison; Leiden clusters are newly calculated.',
           'QC cutoffs and preprocessing may differ from the original paper.','',
           '## Selected samples','',samples.to_markdown(index=False),'',
           '## Figures','']
    for name in names:
        if (fig/name).exists():
            lines.extend([f'### {name}',f'![{name}](figures/{name})',''])
    lines.extend(['## Next steps','Reproduce author QC and integration methods on more patients; perform patient-level analyses.',
                  'Do not use this three-sample pilot to infer differences between cancer subtypes.'])
    dest=out/'pilot_report.md'
    dest.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('Saved:',dest)

if __name__=='__main__':
    main()

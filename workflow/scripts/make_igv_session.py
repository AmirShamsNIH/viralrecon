#!/usr/bin/env python3
"""Write an IGV XML session loading per-sample BAMs, annotated VCFs and consensus
FASTAs beside the reference FASTA for one target."""

import argparse
import os
import xml.etree.ElementTree as ET


def make_session(target, ref, bams, vcfs, fastas, output):
    root = ET.Element(
        "Session",
        attrib={
            "genome":   ref,
            "locus":    "All",
            "version":  "8",
        },
    )

    resources = ET.SubElement(root, "Resources")

    # Reference FASTA
    ET.SubElement(resources, "Resource", path=ref)

    # BAMs
    for bam in bams:
        name = os.path.basename(bam).replace(".bowtie2_map.", ".").split(".bam")[0]
        ET.SubElement(resources, "Resource", path=bam, name=name)

    # VCFs
    for vcf in vcfs:
        name = (
            os.path.basename(vcf)
            .replace(".freebayes_haplotypecaller.", ".")
            .split(".vcf")[0]
        )
        ET.SubElement(resources, "Resource", path=vcf, name=name + " (variants)")

    # Consensus FASTAs
    for fa in (fastas or []):
        name = os.path.basename(fa).split(".consensus.fa")[0]
        ET.SubElement(resources, "Resource", path=fa, name=name + " (consensus)")

    # Panel layout: group tracks by sample (BAM, VCF, then consensus FASTA)
    panel = ET.SubElement(root, "Panel", name="DataPanel")

    for bam in bams:
        sample = os.path.basename(bam).split(".bowtie2_map.")[0]
        ET.SubElement(
            panel,
            "Track",
            attrib={
                "id":          bam,
                "name":        sample,
                "type":        "ALIGNMENT",
                "displayMode": "SQUISHED",
                "visible":     "true",
            },
        )

    for vcf in vcfs:
        sample = os.path.basename(vcf).split(".freebayes_haplotypecaller.")[0]
        ET.SubElement(
            panel,
            "Track",
            attrib={
                "id":          vcf,
                "name":        sample + " variants",
                "type":        "VARIANT",
                "displayMode": "EXPANDED",
                "visible":     "true",
            },
        )

    for fa in (fastas or []):
        sample = os.path.basename(fa).split(".consensus.fa")[0]
        ET.SubElement(
            panel,
            "Track",
            attrib={
                "id":          fa,
                "name":        sample + " consensus",
                "type":        "SEQUENCE",
                "displayMode": "EXPANDED",
                "visible":     "true",
            },
        )

    tree = ET.ElementTree(root)
    ET.indent(tree, space="    ")  # Python 3.9+
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    tree.write(output, xml_declaration=True, encoding="UTF-8")
    print(f"Written: {output}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target",  required=True)
    ap.add_argument("--ref",     required=True)
    ap.add_argument("--bams",    nargs="+", required=True)
    ap.add_argument("--vcfs",    nargs="+", required=True)
    ap.add_argument("--fastas",  nargs="+", default=[],
                    help="Per-sample consensus FASTA files (optional)")
    ap.add_argument("--output",  required=True)
    args = ap.parse_args()
    make_session(args.target, args.ref, args.bams, args.vcfs, args.fastas, args.output)


if __name__ == "__main__":
    main()

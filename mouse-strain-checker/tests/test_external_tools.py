"""Offline alignment smoke test; runs when HISAT2/SAMtools are installed."""
from dataclasses import replace
from pathlib import Path
import random
import shutil
import tempfile
import unittest

from mouse_strain.core import Rules
from mouse_strain.panel import build_panel, load_panel
from mouse_strain.pipeline import build_index, index_identity, map_run, infer


@unittest.skipUnless(all(shutil.which(x) for x in ("hisat2", "hisat2-build", "samtools")),
                     "HISAT2 and SAMtools are required; run inside the supplied Docker image")
class ExternalToolSmokeTest(unittest.TestCase):
    def test_synthetic_fastq_to_hisat2_bam_to_chromosome_candidate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            randomizer = random.Random(82931)
            sequence = "".join(randomizer.choices("ACGT", k=50000))
            chromosome2 = "".join(randomizer.choices("ACGT", k=1000))
            fasta = root / "reference.fa"
            fasta.write_text(">1\n" + sequence + "\n>2\n" + chromosome2 + "\n")
            vcf = root / "strains.vcf"
            fastq = root / "reads.fastq"
            records, reads = [], []
            for pos in range(1000, 49001, 1000):
                ref = sequence[pos - 1]
                alt = next(base for base in "ACGT" if base != ref)
                records.append(f"1\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT:FI\t0/0:1\t1/1:1")
                for number in range(10):
                    start = pos - 40 + number
                    read = list(sequence[start:start + 100])
                    read[pos - 1 - start] = alt
                    reads.append(f"@site_{pos}_fragment_{number}\n{''.join(read)}\n+\n{'I' * 100}\n")
            vcf.write_text("##fileformat=VCFv4.2\n##reference=GRCm39\n"
                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDemo_A\tDemo_B\n"
                           + "\n".join(records) + "\n")
            fastq.write_text("".join(reads))
            panel = root / "panel"
            build_panel(vcf, fasta, panel, "GRCm39")
            prefix = root / "index" / "genome"
            build_index(fasta, prefix, "GRCm39", threads=1)
            meta, db = load_panel(panel)
            db.close()
            idx = index_identity(prefix, meta)
            bam = map_run([fastq], dict(run_accession="SRR0", sample_accession="SRS0"),
                          prefix, idx, root / "mapped.bam", threads=1)
            results = infer(bam, fasta, panel, root / "result",
                            replace(Rules(), window_size=10000))
            self.assertEqual(results[0]["call"], "Demo_B")
            self.assertEqual(results[0]["status"], "candidate_supported")
            self.assertEqual(results[1]["status"], "no_observed_snps")


if __name__ == "__main__":
    unittest.main()
